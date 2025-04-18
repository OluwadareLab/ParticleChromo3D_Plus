from Swarm import Swarm
import Helper

import sys
import argparse
import os
import logging

import numpy as np
from scipy import stats

from multiprocessing import Pool
from multiprocessing import cpu_count
from itertools import repeat

from Helper import Write_Log

# email imports
import smtplib
import uuid
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header

# Initialize logger
logger = logging.getLogger(__name__)

# I should refactor these but I dont think the interest
randRange = None
swarmSize = None
ittCount = None


# Function defs
def lossFunction(tar, b, lossFunctionChoice):
    if lossFunctionChoice == 2:
        newCost = np.sqrt(np.sum((b - tar) ** 2))  # RMSE
    elif lossFunctionChoice == 0:
        newCost = np.sum((b - tar) ** 2)  # SSE
    elif lossFunctionChoice == 1:
        newCost = np.square(np.subtract(b, tar)).mean()
    elif lossFunctionChoice == 3:
        # Heuber
        alpha = 0.5
        y = tar
        yHat = b
        newCost = np.sum(
            np.where(
                np.abs(y - yHat) < alpha,
                0.5 * (y - yHat) ** 2,
                alpha * (np.abs(y - yHat) - 0.5 * alpha),
            )
        )
    else:
        newCost = np.sqrt(np.sum((b - tar) ** 2))  # RMSE

    return newCost


# Prints statistics of the current swarm
def Print_Stats(swarm, contact, pointCount, i, outFilePtr, convFact):
    pers = stats.pearsonr(swarm.gBest[2], contact[:, 3])
    spear = stats.spearmanr(swarm.gBest[2], contact[:, 3])
    spearIF = stats.spearmanr(swarm.gBest[2], contact[:, 2])

    error = np.sqrt((1 / pointCount) * np.sum((swarm.gBest[2] - contact[:, 3]) ** 2))

    logger.debug('id: ' + str(swarm.id) + 
        ' itt: ' + str(i) + 
        ' Cost: ' + str(swarm.gBest[1]) + 
        ' Pearson: ' + str(pers[0]) + 
        ' Spearmen: ' + str(spear[0]) +
        ' IFSpear: ' + str(spearIF[0]) +
        ' error: ' + str(error))


def Write_Stats(swarm, contact, outFilePtr):
    Helper.Write_Output(outFilePtr, swarm.gBest[0])


# Performs one operation and prints statistics of current swarm
def One_Move(
    ittCount,
    swarm,
    contact,
    pointCount,
    threshold,
    outFilePtr,
    convFact,
    lossFunctionChoice,
):
    saveGBestCost = float("inf")
    totTime = 0

    for i in range(ittCount):
        if (i % 1000 == 0) and (swarm.gBest is not None):
            error = lossFunction(contact[:, 3], swarm.gBest[2], lossFunctionChoice)
            Print_Stats(swarm, contact, pointCount, i, outFilePtr, convFact)

            if (np.abs(saveGBestCost - error)) >= threshold:
                saveGBestCost = error
            else:
                return i, totTime

        operation(i, swarm)

    return i


# Performs a single PSO pass: Velocity calculation, update position, get new cost
def operation(i, swarm):
    swarm.Calc_Vel(ittCount, i)
    swarm.Update_Pos(i)
    swarm.Cost()


# Optimizes single swarm
def Optimize(
    inFilePtr,
    outFilePtr,
    convFact,
    constraint,
    points,
    zeroInd,
    randRange,
    swarmSize,
    threshold,
    ittCount,
    lossFunctionChoice,
):
    dist = 1.0 / (constraint[:, 2] ** convFact)
    constraint = np.insert(constraint, 3, dist, axis=1)

    swarm = Swarm(
        constraint, len(points), randVal=randRange, swarmSize=swarmSize, zeroInd=zeroInd
    )

    ittFin = One_Move(
        ittCount,
        swarm,
        constraint,
        len(points),
        threshold,
        outFilePtr,
        convFact,
        lossFunctionChoice,
    )
    return (
        stats.pearsonr(swarm.gBest[2], constraint[:, 3])[0],
        stats.spearmanr(swarm.gBest[2], constraint[:, 3])[0],
        lossFunction(constraint[:, 3], swarm.gBest[2], lossFunctionChoice),
        ittFin,
        swarm.id,
        swarm,
    )


# Runs in paralel if passed multiple rangeSpace
def Par_Choice(inFilePtr, outFilePtr, alpha):
    global randRange, swarmSize, threshold, ittCount, lossFunctionChoice
    contact, points, zeroInd = Helper.Read_Data(inFilePtr, alpha)

    bestSwarm = None
    if 1 == 1:
        convStore = []
        alphas = np.array(range(int(alpha[0]), int(alpha[1]), int(alpha[2]))) / 100
        pool = Pool(processes=PROC_COUNT)
        swarms = pool.starmap(
            Optimize,
            zip(
                repeat(inFilePtr),
                repeat(outFilePtr),
                alphas,
                repeat(contact),
                repeat(points),
                repeat(zeroInd),
                repeat(randRange),
                repeat(swarmSize),
                repeat(threshold),
                repeat(ittCount),
                repeat(lossFunctionChoice),
            ),
        )

        pool.close()
        pool.join()

        iforapl = 0
        for swarm in swarms:
            # print(str(swarm[-1]) + ' ' + str(swarm[1]))
            convStore.append(swarm)
            if (bestSwarm is None) or (swarm[1] > bestSwarm[1]):
                bestSwarm = swarm
                swarmForPDB = swarm[len(swarm) - 1]
                bestAlpha = alphas[iforapl]
            iforapl += 1
    else:
        bestSwarm = Optimize(inFilePtr, outFilePtr, alpha)
    contact = np.insert(contact, 3, 1.0 / (contact[:, 2] ** bestAlpha), axis=1)
    logger.info(f"bestSwarm : {bestSwarm}")
    Write_Stats(swarmForPDB, contact, outFilePtr)

    return bestSwarm


def Full_List(inputFilePtr, outFilePtr, alpha):
    global randRange
    convStore = []

    convStore.append(Par_Choice(inputFilePtr, outFilePtr, alpha))
    logger.info(
        f"pearson: {convStore[0][0]}"
        + f" spearman: {convStore[0][1]}"
        + f" rmse: {convStore[0][2]}"
    )

    # Helper.Write_List(convStore, outFilePtr)
    return convStore

def setup_logging(loglevel):
    level = getattr(logging, loglevel.upper())  # Convert string to corresponding logging level
    logging.basicConfig(level=level,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    if cpu_count() <= 2:
        PROC_COUNT = cpu_count()
    else:
        PROC_COUNT = cpu_count() - 1  # attempt to reduce thrashing.

    rangeSpace = []  # Max scaling factor. Needs to be optimized for each specific dataset. Use two values [one, two] to multithread through a range of those two values at a interval of 5000

    # Arguments for running program
    # python3 ParticleChromo3D.py <input_data> <other_parameter>
    parser = argparse.ArgumentParser("ParticleChromo3D")
    parser.add_argument("infile", help="Matrix of contacts", type=str)
    parser.add_argument(
        "-ss",
        "--swarmSize",
        help="Number of particles in system [Default 15]",
        type=int,
        default=5,
    )
    parser.add_argument(
        "-itt",
        "--ittCount",
        help="Maximum itterations before stop [Default 30000]",
        type=int,
        default=30000,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="Error threshold before stoping [Default 0.000001]",
        type=float,
        default=0.000001,
    )
    parser.add_argument(
        "-rr",
        "--randRange",
        help="Range of x,y,z starting coords. Random value bewtween -randRange,randRange [Default 1]",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Filename of the output pdb model  [Default ./chr.pdb]",
        type=str,
        default="./out/chr.pdb",
    )
    parser.add_argument(
        "-e",
        "--email",
        help="Email to message [Default NULL]",
        type=str,
        default="NULL",
    )
    parser.add_argument(
        "-lf",
        "--lossFunction",
        help="Email to message [Default NULL]",
        type=str,
        default="2",
    )
    parser.add_argument(
        "-ll",
        "--logLevel",
        help="Log level [Default INFO]",
        type=str,
        default="INFO",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],  # Valid log levels
    )

    

    args = parser.parse_args()

    if args.infile:
        inFilePtr = args.infile
    if args.outfile:
        outFilePtr = args.outfile + str(uuid.uuid4())
    else:
        outFilePtr = "noIn"
    if args.swarmSize:
        swarmSize = args.swarmSize
    if args.ittCount:
        ittCount = args.ittCount
    if args.threshold:
        threshold = args.threshold
    randRange = args.randRange
    if args.email:
        emailAddr = args.email
    if args.lossFunction:
        lossFunctionChoice = int(args.lossFunction)

    setup_logging(args.logLevel)  # Set logging level
    logger.info(f"Starting ParticleChromo3D with log level: {args.logLevel}")

    if len(rangeSpace) == 0:
        rangeSpace.append(20000)

    if len(rangeSpace) > 2 and (rangeSpace[0] == rangeSpace[1]):
        rangeSpace.pop()

    logger.info(f"proccessing file : {inFilePtr}")

    fout = inFilePtr + ".stripped"
    clean_lines = []
    f = open(inFilePtr, "r")
    lines = f.readlines()
    for line in lines:
        res = str(" ".join(line.split()))
        clean_lines.append(res)
    f.close()

    with open(fout, "w") as f:
        f.writelines("\n".join(clean_lines))
    f.close()

    theseAlphas = np.array([0.1, 2.0, 0.1]) * 100
    theAlphas = (
        np.array(range(int(theseAlphas[0]), int(theseAlphas[1]), int(theseAlphas[2])))
        / 100
    )

    if outFilePtr == "noIn":
        outFilePtr = os.path.basename(os.path.basename(inFilePtr) + str(uuid.uuid4()))
        outFilePtr = os.path.splitext(outFilePtr)[0]
        logger.info(outFilePtr)
    outputOfSwarm = Full_List(inFilePtr + ".stripped", outFilePtr, theseAlphas)[0]
    logger.info(outputOfSwarm)

    bestSpearm = outputOfSwarm[1]
    bestCost = outputOfSwarm[2]
    bestAlpha = theAlphas[outputOfSwarm[4]]
    bestPearsonRHO = outputOfSwarm[0]

    logger.info(f"Input file: {inFilePtr}")
    logger.info(f"Convert factor:: {bestAlpha}")
    logger.info(f"SSE at best spearman : {bestCost}")
    logger.info(f"Best Spearman correlation Dist vs. Reconstructed Dist  : {bestSpearm}")
    logger.info(f"Best Pearson correlation Dist vs. Reconstructed Dist: {bestPearsonRHO}")
    Write_Log(
        outFilePtr + ".log", inFilePtr, bestAlpha, bestCost, bestSpearm, bestPearsonRHO
    )

    if "HOSTNAME_BE" in os.environ:
        #'<br><br>neo4j by first going to https://biomlearn.uccs.edu:5000/neo and then https://biomlearn.uccs.edu:7474/browser/')
        print(
            "<br><br>Download pdb at: http://"
            + os.environ["HOSTNAME_BE"]
            + ":5001/download?ofname="
            + outFilePtr
        )

    ############################## email section

    if "SERVICE_EMAIL_KEY" not in os.environ or "SERVICE_EMAIL" not in os.environ:
        logger.warning("Missing email properties in PS.py")
        sys.exit(0)

    gmail_pass = os.environ["SERVICE_EMAIL_KEY"]
    user = os.environ["SERVICE_EMAIL"]
    host = "smtp.gmail.com"
    port = 465

    to = emailAddr
    subject = "ParticleChromo3D"
    emailBody = (
        "\n\nAfter processing the file: "
        + inFilePtr
        + " The best Spearman correlation Dist vs. Reconstructed Dist found was : "
        + str(bestSpearm)
        + ". The Best Pearson correlation Dist vs. Reconstructed Dist was : "
        + str(bestPearsonRHO)
        + ".\n\n"
    )
    inputsEBody = (
        "The parameters for this run were set to {ifname="
        + inFilePtr
        + ", ss="
        + str(swarmSize)
        + ", itt="
        + str(ittCount)
        + ", threshold="
        + str(threshold)
        + ", randRange="
        + str(randRange)
        + "} "
    )

    body = (
        "Your ParticleChromo3D+ job completed.\n The results can be found at : "
        + ": http://biomlearn.uccs.edu:5001/download?ofname="
        + outFilePtr
        + emailBody
        + inputsEBody
    )
    filename = "run.sh"

    # locate and attach desired attachments
    fileNameForEmail = outFilePtr + ".pdb"
    f = open(fileNameForEmail, "r")
    att = MIMEText(f.read())
    att.add_header("Content-disposition", "attachment", filename=fileNameForEmail)
    f.close()

    # create message object
    message = MIMEMultipart()

    # add in header
    message["From"] = Header(user)
    message["To"] = Header(to)
    message["Subject"] = Header(subject)

    # attach message body as MIMEText
    message.attach(MIMEText(body, "plain", "utf-8"))

    message.attach(att)

    #  setup email server
    server = smtplib.SMTP_SSL(host, port)
    server.login(user, gmail_pass)

    # send email and quit server
    server.sendmail(user, to, message.as_string())
    server.quit()
