#################################################################
# @Program: utils.py                                            #
# @Version: 1                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 401 Terry Ave North                                           #
# Seattle, Washington  98109-5234                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  5/17/2012                      #
#################################################################

import re

########################################################
# Methods needed to convert miRNA names to miRBase IDs #
########################################################
def miRNAInDict(miRNA, miRNAIDs):
    retMe = []
    for i in miRNAIDs.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        re1 = re.compile(a+'[a-z]$')
        if re1.match(b):
            return 1
    else:
        re1 = re.compile(b+'[a-z]$')
        if re1.match(a):
            return 1
    return 0

def uniquify(array1):
    inThere = []
    outThere = []
    for i in array1:
        if i[0]<i[2]:
            if not i[0]+i[2] in inThere:
                inThere.append(i[0]+i[2])
                outThere.append(i)
        else: 
            if not i[2]+i[0] in inThere:
                inThere.append(i[2]+i[0])
                outThere.append(i)
    return outThere

def randPSSMClustSize(clustSize):
    if clustSize <= 5:
        return 5
    elif clustSize <= 10:
        return 10
    elif clustSize <= 15:
        return 15
    elif clustSize <= 20:
        return 20
    elif clustSize <= 25:
        return 25
    elif clustSize <= 30:
        return 30
    elif clustSize <= 35:
        return 35
    elif clustSize <= 40:
        return 40
    elif clustSize <= 45:
        return 45
    elif clustSize <= 50:
        return 50
    elif clustSize <= 55:
        return 55
    elif clustSize <= 60:
        return 60
    else:
        return 65

# Build a pssms dictionary for TOMTOM
def compilePssms(clusterMemeRuns,maxEValue):
    pssmsNames = []
    pssms = []
    for i in range(len(clusterMemeRuns)):
        motifSet = clusterMemeRuns[i]
        for j in range(len(motifSet)):
            pssm = motifSet[j]
            if float(pssm.getEValue())<=maxEValue:
                pssmsNames.append(str(i)+'.'+str(j))
                pssm.setName(str(i)+'.'+str(j))
                pssms.append(deepcopy(pssm))
    return dict(zip(pssmsNames,pssms))

# Returns a list of all p-values
def getPValues(tomtomResults):
    pValues = []
    allScores = tomtomResults.getAllScores()
    o1 = allScores.keys()[0]
    for i in allScores[o1]:
        pValues.append(allScores[o1][i]['pValue'])
    return pValues


# Make the files for a TomTom run
def makeFiles(nucFreqs, queryPssms, targetPssms, num, strands='+ -'):
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 3.0\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: '+strands+'\n\n'
    memeHeader += 'Background letter frequencies (from genome):\n'
    memeHeader += 'A '+str(round(float(nucFreqs['A']),3))+' C '+str(round(float(nucFreqs['C']),3))+' G '+str(round(float(nucFreqs['G']),3))+' T '+str(round(float(nucFreqs['T']),3))+'\n\n'
    # Make query PSSM file
    queryFile = open('tmp/query'+str(num)+'.tomtom','w')
    queryFile.write(memeHeader)
    queryFile.write('\n\n'.join([pssm1.getMemeFormatted() for pssm1 in queryPssms]))
    queryFile.close()
    # Make target PSSM file
    targetFile = open('tmp/target'+str(num)+'.tomtom','w')
    targetFile.write(memeHeader)
    targetFile.write('\n\n'.join([pssm1.getMemeFormatted() for pssm1 in targetPssms]))
    targetFile.close()

