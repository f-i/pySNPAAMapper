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

inFile = sys.argv[1]  # Stores file one for input checking; note that sys.argv[0] is this script file name

outFile = inFile + ".exons"

#number of bases on exon boundary
exonboundary_offset=1
up_flank=2000
down_flank=500
