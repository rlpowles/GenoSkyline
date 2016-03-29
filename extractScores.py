#!/usr/bin/env python

from __future__ import print_function
from intervaltree import Interval, IntervalTree
from collections import defaultdict
import argparse, progressbar

def getChrNum(Chr):
    """ Sort by chromosome """
    if Chr:
        New = Chr[3:]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New

def calculateWindowLength(positionTuple, chr, windowSize):
    """
    check that the window doesn't go past the beginning or end of chromosome. If it does the denominator of the posterior
    weighted mean should change. Returns the adjusted interval size if the SNP is near the beginning or end of a chromosome

    """
    chrSizes = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
                141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
                81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16571]
    currentChrSize = chrSizes[getChrNum(chr)-1]
    startGap = (abs(positionTuple[0]) + 1) if positionTuple[0] <= 0 else 0
    endGap = (positionTuple[1] - currentChrSize) if positionTuple[1] > currentChrSize else 0
    return(windowSize - startGap - endGap)


description = """

Script for extracting GS Scores from a GenoSkyline annotation bed file at the specific positions of a supplied GWAS dataset.

"""
parser = argparse.ArgumentParser(description=description)
parser.add_argument('gwasPos', help='GWAS positions')
parser.add_argument('annotBed', help='Annotation Bed File')
parser.add_argument('-o', dest="outName", default = "tissueData", help="Output file name")
parser.add_argument('-i', dest='intervalSize', type=int, default = 10000, help="size of interval")
args = parser.parse_args()

#If the window is odd than both left and right sizes are the same, if even, left window is -1
leftwindow = int(args.intervalSize/2 - 1 if (args.intervalSize % 2 ==0) else args.intervalSize/2)
rightwindow = int(args.intervalSize/2)

#Each key is a chromosome, each value is a list of positions in that chromosome
chrDict = defaultdict(list)
print("Reading GWAS data file...")
with open(args.gwasPos) as fp:
    for line in fp:
        c = "chr" + line.split()[0]
        position = int(line.split()[1])
        chrDict[c].append(position)
#Matching intervaltrees containment: start + 1, end + 1
bedTrees = defaultdict(IntervalTree)
with open(args.annotBed) as fp:
    num_lines = sum(1 for line in fp)

print("Loading and indexing bed intervals...")
with open(args.annotBed) as fp:
    for i,line in enumerate(fp):
        c = line.split()[0]
        start = int(line.split()[1])
        stop = int(line.split()[2])
        score = float(line.split()[4])
        bedTrees[c][start+1:stop+1] = score

outfile = open(args.outName, "w")
chrScores = defaultdict(list)
#sort chromosome keys for printing coherency
sortedChrDict = sorted(chrDict.keys(), key=lambda x: getChrNum(x))
for chr in sortedChrDict:
    currentTree = bedTrees[chr]
    print("Calculating scores for {0} positions...".format(chr))
    #get scores for each position with a binary search
    for pos in chrDict[chr]:
        #grab subset of the tree that we'll need for this position
        posStart = pos-leftwindow
        posStop = pos+rightwindow+1
        #Use slice to drop the overhangs of the first and last interval
        currentTree.slice(posStart)
        currentTree.slice(posStop)
        results = currentTree[posStart:posStop]
        normalizedPosterior = [float(i.data)*float(i.end - i.begin) for i in sorted(results)]
        windowLength = calculateWindowLength((posStart, posStop), chr, args.intervalSize)
        chrScores[chr].append(sum(normalizedPosterior)/float(windowLength))
        #check the outside positions to see if they're cutting into intervals
    chrNum = chr.partition("chr")[2]
    for p, s in zip(chrDict[chr], chrScores[chr]):
        print("{0}\t{1}\t{2}".format(chrNum, p, s), file=outfile)
