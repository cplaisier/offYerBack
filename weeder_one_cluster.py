#################################################################
# @Program: weeder_one_cluster.py                               #
# @Version: 1                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  7/5/2010                       #
#################################################################

from optparse import OptionParser
from os import getpid
from copy import deepcopy

#### Option Parsing ####
usage = "%prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-f", "--fasta", dest="fastaFile")
(options, args) = parser.parse_args()

import os
from pssm import pssm
from subprocess import *
from copy import deepcopy
import cPickle

# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
def weeder(seqFile=None, revComp=False):
    if not os.path.exists('tmp/weeder'):
        os.makedirs('tmp/weeder')
    
    # First run weederTFBS
    weederArgs = str(seqFile)+' HS3P small T50'
    if revComp==True:
        weederArgs += ' S'
    #weederArgs += '&> /dev/null'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("progs/weederlauncher " + weederArgs, shell=True,stdout=PIPE, stderr=errOut)
    errOut.close()
    output = weederProc.communicate()

    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqFile)+'.wee','r')
    outLines = [line for line in output.readlines() if line.strip()]
    hitBp = {}
    # Get top hit of 6bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[6] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Searching for motifs of length 8') == -1:
            break

    # Get top hit of 8bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[8] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Your sequences:') == -1:
            break
    
    # Get into the highest ranking motifs
    seqDict = {}
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('**** MY ADVICE ****') == -1:
            break
        splitUp = outLine.strip().split(' ')
        seqDict[splitUp[1]] = splitUp[3].lstrip('>')

    # Get into the highest ranking motifs
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Interesting motifs (highest-ranking)') == -1:
            break
    while 1:
        name = outLines.pop(0).strip() # Get match
        if not name.find('(not highest-ranking)') == -1:
            break
        # Get redundant motifs
        outLines.pop(0)
        redMotifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
        outLines.pop(0)
        outLines.pop(0)
        line = outLines.pop(0)
        instances = []
        while line.find('Frequency Matrix') == -1:
            splitUp = [i for i in line.strip().split(' ') if i]
            instances.append({'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')') })
            line = outLines.pop(0)
        # Read in Frequency Matrix
        outLines.pop(0)
        outLines.pop(0)
        matrix = []
        col = outLines.pop(0)
        while col.find('======') == -1:
            nums = [i for i in col.strip().split('\t')[1].split(' ') if i]
            colSum = 0
            for i in nums:
                colSum += int(i.strip())
            matrix += [[ float(nums[0])/float(colSum), float(nums[1])/float(colSum), float(nums[2])/float(colSum), float(nums[3])/float(colSum)]]
            col = outLines.pop(0)
        PSSMs += [pssm(biclusterName=name,nsites=instances,eValue=hitBp[len(matrix)][1],pssm=matrix,genes=redMotifs)]
    return PSSMs

# Run Weeder
weederPSSMs = weeder(seqFile=options.fastaFile, revComp=False)

# Get the data to make the p-values for the motifs
#pklFile = open('weederRand.pkl','rb')
#weederRand = cPickle.load(pklFile)
#pklFile.close()

# For now just reduce down to 8bp or if no 8bp then a 6bp motif
tmpPSSMs = []
for pssm1 in weederPSSMs:
    if len(pssm1.getName())==8:
        tmpPSSMs.append(deepcopy(pssm1))
if len(tmpPSSMs)==0 and len(weederPSSMs)>0:
    tmpPSSMs.append(weederPSSMs[0])
weederPSSMs = deepcopy(tmpPSSMs)

# Go over each motif and output the files so R can read in the data
motifId = 1
for motif in weederPSSMs:
    width = str(len(motif.getConsensusMotif()))
    llr = 'NA'
    # Calculate e-value
    e_value = motif.getEValue()
    # Get number of targets in genes
    sites = []
    for site in motif.getNSites():
        if not site['gene'] in sites:
            sites.append(site['gene'])
    sites = len(sites)
    outFile1 = open(str(options.fastaFile)+'.'+str(motifId)+'.f1','w')
    outFile1.write(str(width)+','+str(llr)+','+str(e_value)+','+str(sites))
    for i in motif.getMatrix():
        outFile1.write('\n'+str(i[0])+','+str(i[1])+','+str(i[2])+','+str(i[3]))
    outFile1.close()
    outFile2 = open(str(options.fastaFile)+'.'+str(motifId)+'.f2','w')
    cnt = 1
    outFile2.write(',gene,strand,start,p.value,site')
    for site in motif.getNSites():
        outFile2.write('\n'+str(cnt)+','+str(site['gene'])+','+str(site['strand'])+','+str(site['start'])+','+str(site['match'])+','+str(site['site']))
        cnt += 1
    outFile2.close()
    motifId += 1

# Write the MEME formatted file
mastHeader = 'ALPHABET= ACGT\n'
queryFile = open(str(options.fastaFile)+'.meme','w')
queryFile.write(mastHeader)
queryFile.write('\n'.join([pssm.getMastFormatted() for pssm in weederPSSMs]))
queryFile.close()

