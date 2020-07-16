#!/usr/bin/env python
'''
Algorithm for generating annotation information structure for each exon (Python
script to parse exons from ChrAll_knownGene.txt file)

USAGE:
    python Algorithm_generating_annotation_exon.py ChrAll_knownGene.txt
'''

import sys

def argsCheck(numArgs):
    '''Checks if in proper number of arguments are passed; gives instructions on
    proper use'''
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print("\nPlease supply ONLY ONE file to parse.\n\nCorrect Syntax is:\n    python " + sys.argv[0] + " ChrAll_knownGene.txt\n\nThe output file will be the input file name with \".exons\" appended.\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

# Uses the function defined above to check if the number of arguments is correct
argsCheck(2)

#number of bases on exon boundary
exonboundary_offset=1
up_flank=2000
down_flank=500


inFile = sys.argv[1]  # Stores file one for input checking; note that sys.argv[0] is this script file name
try:
    with open(inFile,"r") as in_file:
        exons=in_file.read()

        in_file.close()  # safely closes file
except IOError:  #handling Exceptions
    print("Failed to open " + inFile)
    exit(1)

outFile = inFile + ".exons"
try:
    with open(outFile,"w") as out_file:
        header="bin\tname\tchrom\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tExonNumber\texStart\texEnd\tintronStart\tintronEnd\tUTR5Start\tUTR5End\tUTR3Start\tUTR3End\tupstreamStart\tupstreamEnd\tdownstreamStart\tdownstreamEnd\n"

        out_file.close()
except IOError:
    print("Failed to create " + outFile)
    exit(1)
