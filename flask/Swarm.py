import copy

import numpy as np

from scipy.spatial import distance


''' Swarm Cluster for PSO calculations on HiC data to find optimal distances for xyz points of beads
        This creates a cluster of particles, each particle is represented by a list of xyz coordinates, 
        a point of each particle is the individual coordinate (x,y, or z) of each bead. A different velocity is applied
        to each individual coordinate (x, y, or z).

    Swarms are in a cluster format as this allows python to do vector math on the swarm array which is MANY times more efficient
    than doing individual objects of swarms
    
    Developer: 
    David Vadnais dvadnais@uccs.edu
    DMichael Middleton mmiddlet@uccs.edu 12/4/2020

'''

class Swarm:
    id = 0

    def __init__(self, ref, pointCount, randVal = 0.5, swarmSize = 15, zeroInd=None, swarmComb = None , lossFunc = 0):
        '''
        ref: distance matrix from HiC data \n
        pointCount: Number of beads \n
        randVal: random value to calculate initial x,y,z from \n
        swarmSize: number of particles in cluster \n
        zeroInd: Used to delete in distance calculation \n
        swarmComb: used to combine multiple swarms if doing multiple passes \n
        lossFunctionChoice: chosenLossFunc
        '''

        Swarm.id += 1
        self.id = Swarm.id # ID of swarm if using multiprocessing with multiple swarm clusters

        self.pc = pointCount # pointCount

        # Random Values to get initial xyz coordinates from
        self.randMax = randVal
        self.randMin = -randVal
        self.lossFunc = lossFunc

        self.gBest = None # Global best position

        self.ref = ref # Reference distance matrix
        self.zeroInd = zeroInd

        # Creates a list of particles where each particle is a list of xyz coordinates of size pointCount
        tempList = []
        if swarmComb is None:
            for i in range(swarmSize):
                tempList.append(self.Rand_Cur())
        else:
            for swarm in swarmComb:
                tempList.append(swarm.gBest[0])

            for i in range(swarmSize-len(swarmComb)):
                cutSize = np.random.randint(1, (int)(self.pc/2))
                randCopy = np.random.choice(swarmComb).gBest[0]
                tempPos, mask = self.Rand_shift(randCopy, cutSize, 0.1)
                tempList.append(tempPos)

        self.pos = np.asarray(tempList) # Turns the list into a 3D matrix for vector operations
        self.posBest = copy.copy(self.pos) # Best position of each individual particle
        self.costBest = np.full((self.pos.shape[0],1), np.inf) # Best cost of each individual particle

        self.vel = np.zeros((self.pos.shape)) # Velocity of each particle 
        self.cost = np.full((self.pos.shape[0],1), np.inf) # Cost of each particle

        self.dist = np.zeros((self.pos.shape[0],self.ref.shape[0])) # Distance of each particle

        self.locOpCount = np.zeros((self.pos.shape[0],1)) # Used to check if a particle is not changing after a certain amount of itterations

        self.Cost() # Gets first cost calculations

    def Rand_Cur(self):
        ''' Gets a array of xyz coordinates of size pointCount,3 between a minimum and maximum value '''
        return np.random.uniform(self.randMin,self.randMax, size=(self.pc, 3))

    def Rand_shift(self, copyPos, cutSize, threshold=0.1):
        '''
        Performs a random shift of all the xyz positions in a passed in position array. \n

        Args: 
            copyPos: position array to shift  
            cutSize: Number of values to shift 
            threshold: range of the random shift applied to the positions; values will be uniformly sampled from [-threshold, threshold]
        '''

        temp = copy.copy(copyPos)

        # random boolean mask for which values will be changed
        falseVal = np.zeros(cutSize).astype(bool)
        trueVal = np.ones(temp.shape[0]-cutSize).astype(bool)

        randBool = np.append(falseVal, trueVal)
        np.random.shuffle(randBool)

        mask =  np.tile(randBool.transpose(), (temp.shape[1], 1)).transpose()

        # random matrix the same shape of your data
        r = np.random.uniform(-threshold, threshold, size=temp.shape)

        temp[mask] += r[mask]

        return temp, mask

    def Update_Cost(self, newCost):
        ''' Updates the cost of each particle in cluster '''
        tmpMsk = newCost > self.cost

        self.locOpCount[tmpMsk] += np.ones(self.locOpCount.shape)[tmpMsk]

        self.cost = newCost

        tmpMsk = self.cost < self.costBest
        tmpMsk = tmpMsk.reshape(tmpMsk.shape[0])

        self.posBest[tmpMsk] = copy.copy(self.pos)[tmpMsk]
        self.costBest[tmpMsk] = self.cost[tmpMsk]

        currentBest = np.argmin(self.cost)

        if (self.gBest is None) or (self.cost[currentBest] < self.gBest[1]):
            self.gBest = (copy.copy(self.pos[currentBest]),copy.copy(self.cost[currentBest][0]), copy.copy(self.dist[currentBest]))

    
    def Update_Pos(self, itt):
        '''
        Updates the position of each particle \n
        X(t+1) = X(t)+V(t+1) \n
        
        Args:
            itt: the current iteration of program
        '''

        # Gets a random cutize
        cutSize = np.random.randint(1, (int)(self.pc-1))

        # Adjust threshold depending on iteration to gradually adjust
        if (itt > 500):
            thresh = (1/itt)*100
        else:
            thresh = 1

        # Gets a boolean array of values that checks if a certain number of particles has changed in the last iteration
        # This certain number is calculated in Calc_Const and is based on the number of iterations
        # The boolean array that is calculated is then used to shift any particle that have been calculated to not been changing
        tmpMsk = self.locOpCount > self.Calc_Const(10000,itt,5,15)
        tmpMsk = tmpMsk.reshape(self.locOpCount.shape[0])
        changeInd = np.where(tmpMsk)[0]

        # For each particle that is being changed
        for i in changeInd:
            # If iterations are over a certain amount a full new swarm is calculated
            if(itt < 1000):
                self.pos[i] = self.Rand_Cur()
            else:
                # Otherwise, shift only the cutsize amount of particles			
                self.pos[i], mask = self.Rand_shift(self.pos[i], cutSize, thresh)
                self.vel[i][mask] = np.zeros((self.vel.shape[1], self.vel.shape[2]))[mask]

        # Changes all shifted nodes to have a velocity of zero
        self.vel[tmpMsk] = np.zeros(self.vel.shape)[tmpMsk]
        # Resets local optima counter for changed particle
        self.locOpCount[tmpMsk] = -1

        # Adjusts new positions with velocitys
        self.pos[~tmpMsk] = self.pos[~tmpMsk]+self.vel[~tmpMsk]

    def Calc_Dist(self):
        ''' Calculates the euclidean distance for each particle '''
        for i in range(self.pos.shape[0]):
            if self.zeroInd is None:
                self.dist[i] = distance.pdist(self.pos[i])
            else:
                tmpDist = distance.pdist(self.pos[i])
                tmpDist = np.delete(tmpDist, self.zeroInd)
                self.dist[i] = tmpDist

        return self.dist
    
    def Cost(self):
        ''' Calculates the cost of each particle in the swarm cluster '''
        self.Calc_Dist()
        #print("which lf?" )
        #print(self.lossFunc)
        if (self.lossFunc == 1): 
            newCost =  (1/self.pc) * np.sum( (self.dist-self.ref[:,3])**2, axis=1 ) /len(self.ref[:,3])#MSE
        elif(self.lossFunc == 2): 
            #print("loss func : RMSE" )
            newCost = newCost = np.sqrt(np.sum( (self.dist-self.ref[:,3])**2, axis=1 ))#RMSE
        elif(self.lossFunc == 3): 
            #Heuber
            delta = 0.1
            y = self.ref[:,3]
            yHat = self.dist
            newCost = np.sum(np.where(np.abs(y-yHat) < delta,.5*(y-yHat)**2 , delta*(np.abs(y-yHat)-0.5*delta)),axis=1)
            #print("ME:" ,newCost)
            #print("notME",dnewCost)
        else :
            newCost = np.sum( (self.dist-self.ref[:,3])**2, axis=1 )#SSE
        
        newCost = newCost.reshape((self.pos.shape[0],1))

        self.Update_Cost(newCost)

        return self.cost

    def Calc_Vel(self, ittMax, itt):
        '''
        Calculates the velocitys of each particle \n
        V(t+1) = weight*V(t) + conP*ranP(pBest-X(t)) + conG*ranG(gBest-x(t))
        '''

        ranP = np.random.rand(self.pos.shape[0],self.pos.shape[1],self.pos.shape[2])#local rand  variation
        ranG = np.random.rand(self.pos.shape[0],self.pos.shape[1],self.pos.shape[2])#global rand variarion

        # Way to calculate constants, not used anymore but for future reference
        #weight = self.Calc_Const(ittMax,itt,0.1,0.5)
        #conP = self.Calc_Const(ittMax,itt,0.1,0.5)
        #conG = self.Calc_Const(ittMax,itt,0.1,2)

        # Constants found to work the best
        weight = 0.5
        conP=0.3 #Local confidence
        conG=2.5 #Global confidence

        self.vel = (weight*self.vel) + (conP*ranP*(self.posBest-self.pos)) + (conG*ranG*(self.gBest[0]-self.pos))

        return self.vel 

    # 
    def Calc_Const(self, ittMax, k, W_MIN, W_MAX):
        '''
        Function that calculates a constant value

        Args:
            ittMax: maximum number of itterations
            k: current itterations
            W_MIN: constant min
            W_MAX: constant max
        '''
        if ittMax < k:
            k = ittMax
        return ( ( (W_MAX-W_MIN)*((ittMax-k)/ittMax) ) + W_MIN )