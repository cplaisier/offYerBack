#################################################################
# @Program: memeBackgroundDistribution.py                       #
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
# Copyrighted by Chris Plaisier  12/4/2009                      #
#################################################################

import sys, re, os, math, shutil
from pssm import pssm
from tomtom import tomtom
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
from copy import deepcopy
from random import sample
import cPickle

# Run meme and get the output into PSSMs
def meme(num, seqFile=None, bgFile=None, nMotifs=1, minMotifWidth=6, maxMotifWidth=12, revComp=True, seed=None):
    if not os.path.exists('tmp/meme'):
        os.makedirs('tmp/meme')
    # Arguments for tomtom
    memeArgs = str(seqFile)+' -bfile tmp/meme/bgFile.meme -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw ' + str(minMotifWidth) + ' -maxw ' + str(maxMotifWidth) + ' -nmotifs ' + str(nMotifs)
    if revComp==True:
        memeArgs += ' -revcomp'
    if not seed==None:
        memeArgs += ' -cons ' + str(seed)
    print memeArgs
    #errOut = open('tmp/meme/stderr.out','w')
    memeProc = Popen("meme " + memeArgs, shell=True,stdout=PIPE) #,stderr=errOut)
    output = memeProc.communicate()[0].split('\n')
    #errOut.close()
    
    PSSMs = []
    # Now iterate through output and save data
    for i in range(len(output)):
        splitUp1 = output[i].strip().split(' ')
        if splitUp1[0]=='Motif' and splitUp1[2]=='position-specific' and splitUp1[3]=='probability':
            i += 2 # Skip the separator line, go to the summary line
            splitUp = output[i].strip().split(' ')
            width = int(splitUp[5])
            sites = splitUp[7]
            eValue = splitUp[9]
            matrix = []
            for j in range(width):
                i += 1
                matrix += [[float(let) for let in output[i].strip().split(' ') if let]]
            PSSMs.append(pssm(biclusterName=str(splitUp1[1]),nsites=sites,eValue=eValue,pssm=matrix,genes=[]))
    clusterMemeRuns[num] = PSSMs

# Function to run the meme function
def runMeme(i):
    #print clusterFileNames[i]
    meme(i,seqFile=clusterFileNames[i],bgFile=allVars['bgFile'],nMotifs=allVars['nMotifs'],minMotifWidth=allVars['minMotifWidth'], maxMotifWidth=allVars['maxMotifWidth'], revComp=allVars['revComp'], seed=seeds[i])

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

# Parameters for a run
permutations = 1000 # Number of times to run meme, get nMotifs per permutation
maxEValue = 10
clusterSizes = range(4,65)
#pValueThreshold = 0.05 # Threshold for TOMTOM similarity

# TODO!!! These parameters could be gotten from the cMonkey Run
nMotifs = 2
regions = ['upstream', '3pUTR']
motifWidth = { 'upstream': [6, 12], '3pUTR': [4, 9] }
revComp = { 'upstream': True, '3pUTR': None }

# Load up the cMonkey run to get permuted pValues
from cMonkeyWrapper import cMonkeyWrapper
c1 = cMonkeyWrapper('cmonkey-run-hsa.RData')

# Getting 
allSeqs = c1.getSeqsUpstream()
bgFile = 'tmp/meme/bgFile.meme'
seed = None
print 'Not using a seed.'
        
# Get Sequences for the run
if not allSeqs==None:
    allSeqs = allSeqs
else:
    allSeqs = {}
    asf = open(allSeqsFile,'r')
    for line in asf.readlines():
        splitUp = line.strip().split(',')
        if splitUp[0] in seqsInAnalysis:
            allSeqs[splitUp[0]] = splitUp[1]

# Make background distribution for meme runs
if not os.path.exists(bgFile):
    print 'Making bacground distribution...'
    getBackground(allSeqs.values())
    print 'Done.\n'

# Iterate through the regions
pssmDepot = { 'upstream': {}, '3pUTR': {} }
for region in regions:
    # Iterate through each cluster size
    for size in clusterSizes:
        if not os.path.exists('tmp/meme/fasta'):
            os.makedirs('tmp/meme/fasta')
        print region, size
        #seed = orig.getConsensusMotif()

        # Make clusterMembership which is a list of lists of genes
        geneNames = allSeqs.keys()
        clusterMembership = []
        for i in range(permutations):
            clusterMembership.append(sample(geneNames,size))

        # Prepare for running meme for each cluster in clusterMembership
        print 'Starting runs for '+str(len(clusterMembership))+' clusters...'
        mgr = Manager()
        clusterFileNames = mgr.list()
        clusterMemeRuns = mgr.list([i for i in range(len(clusterMembership))])
        if not seed==None:
            seeds = mgr.list(seed)
        else:
            seeds = mgr.list([None for i in range(len(clusterMembership))])
        allVars = mgr.dict( { 'bgFile': bgFile, 'nMotifs': nMotifs, 'minMotifWidth': motifWidth[region][0], 'maxMotifWidth': motifWidth[region][1], 'revComp': revComp[region] } )
        for i in range(len(clusterMembership)):
            cluster = clusterMembership[i]
            clusterFileName = 'tmp/meme/fasta/cluster'+str(i)+'.fasta'
            clusterFileNames.append(clusterFileName)
            memeFile = open(clusterFileName,'w')
            seqs = []
            for gene in cluster:
                if gene in allSeqs:
                    if not allSeqs[gene] in seqs:
                        seqs.append(allSeqs[gene])
                        memeFile.write('>'+str(gene)+'\n'+str(allSeqs[gene])+'\n')
            memeFile.close()

        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs avialable.' 
        pool = Pool(processes=cpus)
        pool.map(runMeme,[i for i in range(len(clusterMembership))])
        print 'Done with clusters.\n'

        # Save the pssms
        similar = 0
        pssms = compilePssms(clusterMemeRuns,maxEValue)
        pssmDepot[region][size] = deepcopy(pssms)
        output = open('pssms_'+region+'_'+str(size)+'.pkl','wb')
        cPickle.dump(pssms,output)
        output.close()

        # Clean-up and prepare for next run
        shutil.rmtree('tmp/meme/fasta')

# Pickle the data so we can use another script to analyze it
output = open('pssms.pkl','wb')
cPickle.dump(pssmDepot,output)
output.close()

