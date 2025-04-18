import numpy as np
import csv
import sys

def squarify(mat):
    
    maxBin = 0
    for row in mat:
        maxBin = max(maxBin, row[0]+1, row[1]+1) 
        
    squareMat = np.zeros((maxBin,maxBin))
    
    for row in mat:
        squareMat[row[0]][row[1]] = row[2]
        squareMat[row[1]][row[0]] = row[2]
    
    return squareMat

fileloc = sys.argv[1]
threeMat = []
with open(fileloc, 'r') as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        print(row)
        newRow = [int(row[0]), int(row[1]), float(row[2])]
        threeMat.append(newRow)


squareMat = squarify(threeMat)

with open('convert.out', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    for row in squareMat:

        writer.writerow(row)

exit(0)
