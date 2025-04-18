import numpy as np
import csv
import sys
import argparse
from ParticleChromo3D.particle_chromo_logger import setup_logger

def squarify(mat):
    '''
    takes
    [name, bin1, bin2]

    nxn where n is number of bins
    '''
    
    maxBin = 0
    for row in mat:
        maxBin = max(maxBin, row[0]+1, row[1]+1) 
        
    squareMat = np.ones((maxBin,maxBin))
    
    for row in mat:
        squareMat[row[0]][row[1]] = row[2]
        squareMat[row[1]][row[0]] = row[2]
    
    return squareMat


def main(input_file, output, delimiter = '\t'):   
    logger = setup_logger()

    threeMat = []
    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for index, row in enumerate(reader):
            if len(row) != 3:
                logger.warning(f"skipped line {index} with content : {row}")
                continue

            row[0] = row[0].strip() # remove spaces

            logger.debug(f"row output : {row}")
            newRow = [
                int(row[0]), 
                int(row[1]), 
                float(row[2])
            ]
            threeMat.append(newRow)


    squareMat = squarify(threeMat)

    with open(output, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for row in squareMat:
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform VCM file into a square matrix.")
    parser.add_argument("input_file", help="Path to the input file (TSV format)")
    parser.add_argument("-o", "--output", default="convert.out", help="Output file name (default: convert.out)")
    parser.add_argument("-d", "--delimiter", default="\t", help="Delimiter used in the input file (default: tab)")
    
    args = parser.parse_args()

    main(args.input_file,  output=args.output, delimiter=args.delimiter)
    sys.exit(0)