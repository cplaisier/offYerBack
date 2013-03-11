#################################################################
# @Program: offYerBackV2.py                                     #
# @Version: 2                                                   #
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
# Copyrighted by Chris Plaisier  5/21/2012                      #
#################################################################


#################################################################
## Parameters                                                  ##
#################################################################

# For MEME analysis
bgFile = 'seqs/bgFile.meme'
nMotifs = 2
regions = [ 'upstream' ]
motifWidth = { 'upstream': [6, 12] }
revComp = { 'upstream': True }
# Filtering parameters
maxScore = 0
maxResidual = 0.5
maxEValue = 10
maxSurv = 0.05/720
#pValueThreshold = 0.05/((375**2)/2)
postProcessedFile = 'postProcessed_gbmTCGA.csv'
randPssmsDir = 'randPSSMs'


#################################################################
## Functions                                                   ##
#################################################################

# Run meme and get the output into PSSMs
def meme(num, seqFile=None, bgFile=None, nMotifs=1, minMotifWidth=6, maxMotifWidth=12, revComp=True, seed=None):
    # Arguments for meme
    memeArgs = str(seqFile)+' -bfile '+ str(bgFile) +' -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw ' + str(minMotifWidth) + ' -maxw ' + str(maxMotifWidth) + ' -nmotifs ' + str(nMotifs)
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
            PSSMs.append(pssm(biclusterName=str((seqFile.split('_')[1]).split('.')[0])+'_motif'+str(splitUp1[1])+'_meme',nsites=sites,eValue=eValue,pssm=matrix,genes=[], de_novo_method='meme'))
    clusterMemeRuns[num] = PSSMs

# Wrapper function to run the meme function using a multiprocessing pool
def runMeme(i):
    meme(i,seqFile=clusterFileNames[i],bgFile=memeVars['bgFile'],nMotifs=memeVars['nMotifs'],minMotifWidth=memeVars['minMotifWidth'], maxMotifWidth=memeVars['maxMotifWidth'], revComp=memeVars['revComp'])


# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
def weeder(bicluster, seqFile=None, bgFile='HS', size='small', enriched='T50', revComp=False):
    print seqFile

    # First run weederTFBS
    weederArgs = str(seqFile)+' '+str(bgFile)+' '+str(size)+' '+str(enriched)
    if revComp==True:
        weederArgs += ' S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weederlauncher " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
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

    # Scroll to where the 8bp reads will be
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

    if size=='medium':
        # Scroll to where the 10bp reads wll be
        while 1:
            outLine = outLines.pop(0)
            if not outLine.find('Searching for motifs of length 10') == -1:
                break

        # Get top hit of 10bp look for "1)"
        while 1:
            outLine = outLines.pop(0)
            if not outLine.find('1) ') == -1:
                break
        hitBp[10] = outLine.strip().split(' ')[1:]

    # Scroll to where the 10bp reads will be
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
    motifId = 1
    biclusterId = str((seqFile.split('_')[1]).split('.')[0])
    while 1:
        if len(outLines)<=1:
            break
        if revComp==True:
            name = outLines.pop(0).strip() # Get match
            name += '_'+outLines.pop(0).strip()
        if revComp==False:
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
        PSSMs += [pssm(biclusterName=str(biclusterId)+'_motif'+str(motifId)+'_weeder',nsites=instances,eValue=hitBp[len(matrix)][1],pssm=matrix,genes=redMotifs, de_novo_method='weeder')]
        motifId += 1
    weederResults[bicluster] = PSSMs

# Wrapper function to run weeder using a multiprocessing pool
def runWeeder(i):
    weeder(i,seqFile=clusterFileNames[i],bgFile=weederVars['bgFile'], size=weederVars['size'], enriched=weederVars['enriched'], revComp=weederVars['revComp'])


# Run TomTom on the files
def TomTom(num, distMeth='ed', qThresh='1', minOverlap=6):
    # Arguments for tomtom
    tomtomArgs = ' -query tmp/query'+str(num)+'.tomtom -target tmp/target'+str(num)+'.tomtom -dist '+str(distMeth)+' -o tmp/tomtom_out -text -q-thresh '+str(qThresh)+' -min-overlap '+str(minOverlap)+' -verbosity 0'
    print tomtomArgs
    #p = Popen("tomtom" + tomtomArgs, shell=True)
    #sts = os.waitpid(p.pid, 0)
    errOut = open('tmp/stderr_'+str(num)+'.out','w')
    tomtomProc = Popen("tomtom" + tomtomArgs, shell=True,stdout=PIPE, stderr=errOut)
    outputFile = open('tmp/tomtom_out/tomtom'+str(num)+'.out','w')
    output = tomtomProc.communicate()[0]
    outputFile.write(output)
    outputFile.close()
    errOut.close()

# Wrapper function to run TomTom using multiprocessing pool
def runTomTom(i):
    TomTom(i, distMeth='ed', qThresh='1', minOverlap=6) #blic5

def phyper(q, m, n, k):
    # Get an array of values to run
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    for i in range(len(q)):
        runMe.append('phyper('+str(q[i])+','+str(m[i])+','+str(n[i])+','+str(k[i])+',lower.tail=F)')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    return [line.strip().split(' ')[1] for line in out[0].strip().split('\n') if line]

# Get first principle component
def firstPrincipalComponent(matrix):
    """
    Cacluate the first prinicipal component of a gene expression matrix.
    Input: Expression matrix gene (rows) x conditions (columns), expects that the
        python equivalent will be an array of arrays of expression values (rows).
    Returns: Array of first prinicipal component and variance explained by first
        principal component.
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    rbind = 'm1 = rbind('
    rbind += ','.join(['c('+','.join([str(i) for i in row])+')' for row in matrix])
    rbind += ')'
    runMe.append(rbind)
    # Compute first principal component
    runMe.append('tmp.pc = try(princomp(t(m1)),TRUE)')
    runMe.append('pc.1 = NA')
    runMe.append('var.exp = NA')
    runMe.append('if(!class(tmp.pc)==\'try-error\') {')
    runMe.append('    pc.1 = tmp.pc$scores[,1]')
    runMe.append('    var.exp = ((tmp.pc$sdev^2)/sum(tmp.pc$sdev^2))[1]')
    runMe.append('}')
    runMe.append('var.exp')
    runMe.append('pc.1')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    # Process output
    splitUp = out[0].strip().split('\n')
    splitUp.pop(0) # Get rid of header dealie
    varExp = float(splitUp.pop(0).strip())
    pc1 = []
    for r1 in splitUp:
        pc1 += [float(i) for i in r1.split(' ') if i and (not i.count('[')==1)]
    # Return output
    return [pc1, varExp]

# Get a correlation p-value from R
def correlation(a1, a2):
    """
    Calculate the correlation coefficient and p-value between two variables.
    Input: Two arrays of float or integers.
    Returns: Corrleation coefficient and p-value.
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    runMe.append('c1 = cor.test(c('+','.join([str(i) for i in a1])+'),c('+','.join([str(i) for i in a2])+'))')
    runMe.append('c1$estimate')
    runMe.append('c1$p.value')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    # Process output
    splitUp = out[0].strip().split('\n')
    rho = float(splitUp[1])
    pValue = float((splitUp[2].split(' '))[1])
    return [rho, pValue]

# Compute survival p-value from R
def survival(survival, dead, pc1, age):
    """
    Calculate the survival correlation coefficient and p-value between two variables.
    Input: Four arrays of float or integers.
    Returns:
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    runMe.append('library(survival)')
    runMe.append('s1 = c('+','.join([str(i) for i in survival])+')')
    runMe.append('dead1 = c('+','.join(['\''+str(i)+'\'' for i in dead])+')')
    runMe.append('pc1 = c('+','.join([str(i) for i in pc1])+')')
    runMe.append('age1 = c('+','.join([str(i) for i in age])+')')
    runMe.append('scph1 = summary(coxph(Surv(s1,dead1==\'DEAD\') ~ pc1))')
    runMe.append('scph2 = summary(coxph(Surv(s1,dead1==\'DEAD\') ~ pc1 + age1))')
    runMe.append('scph1$coef[1,4]')
    runMe.append('scph1$coef[1,5]')
    runMe.append('scph2$coef[1,4]')
    runMe.append('scph2$coef[1,5]')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    # Process output
    splitUp = out[0].strip().split('\n')
    z1 = float((splitUp[0].split(' '))[1])
    pValue1 = float((splitUp[1].split(' '))[1])
    z2 = float((splitUp[2].split(' '))[1])
    pValue2 = float((splitUp[3].split(' '))[1])
    return [[z1, pValue1], [z2, pValue2]]

# To test the survival function
#survival([10,20,30,10,15], ['DEAD','ALIVE','ALIVE','DEAD','ALIVE'], [0.1,0.25,0.4,0.6,0.9], [25,35,45,55,65])

#################################################################
## Python Modules Loading                                      ##
#################################################################

# Default python libraries
from math import log10
import cPickle, os, re, sys
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
from shutil import rmtree
from copy import deepcopy
# Custom offYerBack libraries
from cMonkeyWrapper import cMonkeyWrapper
from pssm import pssm
from miRvestigator import miRvestigator
from tomtom import tomtom
# Functions needed to run this script
from utils import *
from sys import stdout, exit
#from weeder import *
import gzip


#################################################################
## Preparing directories for output                            ##
#################################################################
if not os.path.exists('output'):
    os.makedirs('output')


#######################################################################
## Create a dictionary to convert the miRNAs to there respective ids ##
#######################################################################
inFile = open('miRNA/hsa.mature.fa','r')  # Need to update this to the latest database
miRNAIDs = {}
miRNAIDs_rev = {}
clusterFileNames = {}
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.split(' ')
    if not splitUp[1] in miRNAIDs_rev:
        miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()
    if not splitUp[0].lower() in miRNAIDs:
        miRNAIDs[splitUp[0].lower()] = splitUp[1]
    else:
        print 'Uh oh!',splitUp

if not os.path.exists('output/c1_all.pkl'):
    #################################################################
    ## Load cMonkey Object - turns cMonkey data into objects       ##
    #################################################################

    # If this is the first time then load from the RData file
    if not os.path.exists('output/c1.pkl'):
        c1 = cMonkeyWrapper('../cmonkey_run.db',meme_upstream=0,weeder_upstream=1,weeder_3pUTR=0,pita_3pUTR=0,targetscan_3pUTR=0)
        pklFile = open('output/c1.pkl','wb')
        cPickle.dump(c1,pklFile)
    # Otherwise load the dumped pickle file if it exists
    else:
        print 'Loading precached cMonkey object..'
        pklFile = open('output/c1.pkl','rb')
        c1 = cPickle.load(pklFile)
        print 'Done.\n'
    pklFile.close()


    #c1.meme_upstream = 1 # Works
    #c1.weeder_upstream = 1 # Works
    #c1.weeder_3pUTR = 1 # Works
    #c1.pita_3pUTR = 1 # Works
    #c1.targetscan_3pUTR = 1 # Works

    #################################################################
    ## Fill in the missing parts                                   ##
    #################################################################
     #  Check to see if all parts are there:                       #
     #   A. Upstream motifs (MEME)                                 #
     #   B. Upstream motif (Weeder)                                #
     #   C. 3' UTR Weeder-miRvestigator (Weeder)                   #
     #   D. 3' UTR PITA (Set Enrichment)                           #
     #   E. 3' UTR TargetScan (Set Enrichment)                     #
     ###############################################################

    ## A. Upstream motifs (MEME) ##
    # If MEME hasn't been run on the biclusters upstream sequences then do so
    if not c1.meme_upstream==1:
        print 'Running MEME on Upstreams:'
        # Use already run MEME results if they exist
        if not os.path.exists('output/meme_upstream.pkl'):
            # Make needed directories
            if os.path.exists('tmp'):
                rmtree('tmp')
            if not os.path.exists('tmp/meme/fasta'):
                os.makedirs('tmp/meme/fasta')

            # Run MEME on all biclusters
            mgr = Manager()
            clusterFileNames = mgr.dict()
            o1 = []

            # First make fasta files for all biclusters
            print 'Making Files for MEME Upstream run...'
            for b1 in c1.getBiclusters():
                clusterFileName = 'tmp/meme/fasta/bicluster_'+str(b1)+'.fasta'
                seqs = c1.getBiclusterSeqsUpstream(b1)
                if len(seqs)>0:
                    o1.append(b1)
                    clusterFileNames[b1] = clusterFileName
                    fastaFile = open(clusterFileName,'w')
                    fastaFile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))
                    fastaFile.close()

            # Where all the results will be stored
            clusterMemeRuns = mgr.dict()

            # Parameters to use for running MEME
            memeVars = mgr.dict( { 'bgFile': bgFile, 'nMotifs': nMotifs, 'minMotifWidth': motifWidth['upstream'][0], 'maxMotifWidth': motifWidth['upstream'][1], 'revComp': revComp['upstream'] } )

            # Then run MEME using all cores available
            print 'Running MEME on Upstream sequences...'
            cpus = cpu_count()
            print 'There are', cpus,'CPUs avialable.'
            pool = Pool(processes=cpus)
            pool.map(runMeme,[i for i in o1])
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            pklFile = open('output/meme_upstream.pkl','wb')
            cPickle.dump(deepcopy(clusterMemeRuns),pklFile)
        else:
            print 'Loading from precached object...'
            pklFile = open('output/meme_upstream.pkl','rb')
            clusterMemeRuns = cPickle.load(pklFile)
        pklFile.close()

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i in clusterMemeRuns.keys():
            for pssm1 in clusterMemeRuns[i]:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('meme')
                b1.addPssmUpstream(pssm1)

        print 'Done with MEMEing.\n'

        # MEME upstream has been run on cMonkey run
        c1.meme_upstream = 1

    ## B. Upstream motifs (Weeder) ##
    # If Weeder hasn't been run on the biclusters upstream sequences then do so
    if not c1.weeder_upstream==1:
        print 'Running Weeder on Upstreams:'
        # If this has been run previously just load it up
        if not os.path.exists('output/weeder_upstream.pkl'):
            # Make needed directories
            if os.path.exists('tmp'):
                rmtree('tmp')
            if not os.path.exists('tmp/weeder/fasta'):
                os.makedirs('tmp/weeder/fasta')

            # Run Weeder on all biclusters
            mgr = Manager()
            clusterFileNames = mgr.dict()
            o1 = []

            # First make fasta files for all biclusters
            print 'Making Files for Upstream Weeder re-run...'
            for b1 in c1.getBiclusters():
                clusterFileName = 'tmp/weeder/fasta/bicluster_'+str(b1)+'.fasta'
                seqs = c1.getBiclusterSeqsUpstream(b1)
                if len(seqs)>0:
                    o1.append(b1)
                    clusterFileNames[b1] = clusterFileName
                    fastaFile = open(clusterFileName,'w')
                    fastaFile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))
                    fastaFile.close()

            # Where all the results will be stored
            weederResults = mgr.dict()

            # Parameters to use for running Weeder
            # Set to run Weeder on 'medium' setting which means 6bp, 8bp and 10bp motifs
            #weederVars = mgr.dict( { 'bgFile': 'HS', 'size': 'medium', 'enriched': 'T50', 'revComp': True } )
            weederVars = mgr.dict( { 'bgFile': 'HS', 'size': 'small', 'enriched': 'T50', 'revComp': True } )

            # Run this using all cores available
            print 'Running Weeder...'
            cpus = cpu_count()
            print 'There are', cpus,'CPUs avialable.'
            pool = Pool(processes=cpus)
            pool.map(runWeeder,[i for i in o1])
            pool.close()
            pool.join()

            #runWeeder(40)
            #for i in o1:
            #    runWeeder(i)
            # Dump weeder results as a pickle file
            pklFile = open('output/weeder_upstream.pkl','wb')
            cPickle.dump(deepcopy(weederResults),pklFile)
        else:
            print 'Loading from precached object...'
            pklFile = open('output/weeder_upstream.pkl','rb')
            weederResults = cPickle.load(pklFile)
        pklFile.close()

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i in weederResults.keys():
            for pssm1 in weederResults[i]:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('weeder')
                b1.addPssmUpstream(pssm1)
        print 'Done with Weedering.\n'

        # Weeder upstream has been run on cMonkey run
        c1.weeder_upstream = 1

    # C. 3' UTR Weeder-miRvestigator (Weeder)
    # If Weeder hasn't been run on the biclusters 3' UTR sequences then do so
    if not c1.weeder_3pUTR==1:
        print 'Running Weeder on 3\' UTRs:'
        # If this has been run previously just load it up
        if not os.path.exists('output/weeder_3pUTR.pkl'):
            # Make needed directories
            if os.path.exists('tmp'):
                rmtree('tmp')
            if not os.path.exists('tmp/weeder/fasta'):
                os.makedirs('tmp/weeder/fasta')

            # Run MEME on all biclusters
            mgr = Manager()
            clusterFileNames = mgr.dict()
            o1 = []

            # First make fasta files for all biclusters
            print 'Making Files for Upstream Weeder re-run...'
            for b1 in c1.getBiclusters():
                clusterFileName = 'tmp/weeder/fasta/bicluster_'+str(b1)+'.fasta'
                seqs = c1.getBiclusterSeqs3pUTR(b1)
                if len(seqs)>0:
                    o1.append(b1)
                    clusterFileNames[b1] = clusterFileName
                    fastaFile = open(clusterFileName,'w')
                    fastaFile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))
                    fastaFile.close()

            # Where all the results will be stored
            weederResults = mgr.dict()

            # Parameters to use for running Weeder
            # Set to run Weeder on 'medium' setting which means 6bp, 8bp and 10bp motifs
            weederVars = mgr.dict( { 'bgFile': 'HS3P', 'size': 'small', 'enriched': 'T50', 'revComp': False } )

            # Run this using all cores available
            print 'Running Weeder...'
            cpus = cpu_count()
            print 'There are', cpus,'CPUs avialable.'
            pool = Pool(processes=cpus)
            pool.map(runWeeder,[i for i in o1])
            pool.close()
            pool.join()
            # Dump weeder results as a pickle file
            pklFile = open('output/weeder_3pUTR.pkl','wb')
            cPickle.dump(deepcopy(weederResults),pklFile)
        else:
            print 'Loading from precached object...'
            pklFile = open('output/weeder_3pUTR.pkl','rb')
            weederResults = cPickle.load(pklFile)
        pklFile.close()

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i in weederResults.keys():
            for pssm1 in weederResults[i]:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('weeder')
                b1.addPssm3pUTR(pssm1)
        print 'Done with 3\'UTR Weedering.\n'

        # Weeder 3pUTR has been run on cMonkey run
        c1.weeder_3pUTR = 1

    # D. 3' UTR PITA
    # If PITA enrichment hasn't been calculated for the biclusters then do so
    if not c1.pita_3pUTR==1:
        print 'Running PITA on Biclusters:'
        # If this has been run previously just load it up
        if not os.path.exists('output/pita_3pUTR.pkl'):
            # Get ready for multiprocessor goodness
            mgr = Manager()

            # Get a list of all genes in the biclusters
            print 'Get a list of all genes in run...'
            tmpDict = c1.getBiclusters()
            genesInBiclusters = []
            for bicluster in tmpDict:
                genes = tmpDict[bicluster].getGenes()
                for gene in genes:
                    if not gene in genesInBiclusters:
                        genesInBiclusters.append(gene)
            biclusters = mgr.dict(tmpDict)
            del tmpDict

            # Load up PITA miRNA ids
            # If this is the first time then load from the RData file
            if not os.path.exists('miRNA/pita.pkl'):
                print 'Loading PITA predictions...'
                tmpList = []
                tmpDict = {}
                inFile = open('miRNA/pita_miRNA_sets_geneSymbol.csv','r')
                inLines = [i.strip().split(',') for i in inFile.readlines() if i.strip()]
                for line in inLines:
                    if line[1].upper() in genesInBiclusters:
                        if not line[1] in tmpList:
                            tmpList.append(line[1].upper())
                        if not line[0] in tmpDict:
                            tmpDict[line[0]] = [line[1].upper()]
                        else:
                            tmpDict[line[0]].append(line[1].upper())
                inFile.close()
                pklFile = open('miRNA/pita.pkl','wb')
                cPickle.dump(tmpDict,pklFile)
                cPickle.dump(tmpList,pklFile)
            # Otherwise load the dumped pickle file if it exists
            else:
                print 'Loading precached PITA predictions...'
                pklFile = open('miRNA/pita.pkl','rb')
                tmpDict = cPickle.load(pklFile)
                tmpList = cPickle.load(pklFile)
            pklFile.close()
            # Setup for analysis
            predDict = mgr.dict(tmpDict)
            pred_totalTargets = mgr.list()
            pred_totalTargets.append(set(tmpList))
            del tmpDict
            del tmpList
            print 'PITA has', len(predDict.keys()),'miRNAs.'

            def clusterHypergeo_pita(biclustId, db = predDict, allGenes = pred_totalTargets[0]):
                # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
                # Take gene list and compute overlap with each miRNA
                genes = allGenes.intersection(biclusters[biclustId].getGenes())
                keys1 = db.keys()
                m1s = []
                q = []
                m = []
                n = []
                k = []
                for m1 in keys1:
                    m1s.append(m1)
                    miRNAGenes = allGenes.intersection(db[m1])
                    q.append(len(set(miRNAGenes).intersection(genes)))
                    m.append(len(miRNAGenes))
                    n.append(len(allGenes)-len(miRNAGenes))
                    k.append(len(genes))
                results = phyper(q,m,n,k)
                min_miRNA = []
                perc_targets = []
                min_pValue = float(1)
                for i in range(1,len(results)):
                    if float(results[i]) <= float(0.05)/float(674) and not q[i]==0 and float(q[i])/float(k[i]) >= 0.1:
                        if min_miRNA==[] or float(results[i]) < min_pValue:
                            min_miRNA = [i]
                            perc_targets = [float(q[i])/float(k[i])]
                            min_pValue = float(results[i])
                        elif float(results[i])==min_pValue:
                            min_miRNA.append(i)
                            perc_targets.append(float(q[i])/float(k[i]))
                print 'Bicluster #', biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA])
                return [biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA]), ' '.join([str(j) for j in perc_targets]), min_pValue]

            # Set allGenes as intersect of pita_totalTargets and genesInBiclusters
            allGenes = pred_totalTargets[0].intersection(genesInBiclusters)

            # Run this using all cores available
            print 'Running PITA enrichment analyses...'
            cpus = cpu_count()
            biclustIds = biclusters.keys()
            pool = Pool(processes=cpus)
            res1 = pool.map(clusterHypergeo_pita,biclustIds)
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            pklFile = open('output/pita_3pUTR.pkl','wb')
            cPickle.dump(res1,pklFile)
        else:
            print 'Loading precached analysis...'
            pklFile = open('output/pita_3pUTR.pkl','rb')
            res1 = cPickle.load(pklFile)
        pklFile.close()

        # Stuff into biclusters
        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += miRNAInDict(m1.lower(),miRNAIDs)
            b1.addAttribute('pita_3pUTR',{'miRNA':r1[1], 'percentTargets':r1[2], 'pValue':r1[3], 'mature_seq_ids':miRNA_mature_seq_ids })
        print 'Done.\n'

        # PITA 3pUTR has been run on cMonkey run
        c1.pita_3pUTR = 1

    # E. 3' UTR TargetScan
    # If TargetScan enrichment hasn't been calculated for the biclusters then do so
    if not c1.targetscan_3pUTR==1:
        print 'Running TargetScan on Biclusters:'
        # If this has been run previously just load it up
        if not os.path.exists('output/targetscan_3pUTR.pkl'):
            # Get ready for multiprocessor goodness
            mgr = Manager()

            # Get a list of all genes in the biclusters
            print 'Get a list of all genes in run...'
            tmpDict = c1.getBiclusters()
            genesInBiclusters = []
            for bicluster in tmpDict:
                genes = tmpDict[bicluster].getGenes()
                for gene in genes:
                    if not gene in genesInBiclusters:
                        genesInBiclusters.append(gene)
            biclusters = mgr.dict(tmpDict)
            del tmpDict

            # Load up TargetScan miRNA ids
            if not os.path.exists('miRNA/targetScan.pkl'):
                print 'Loading TargetScan predictions...'
                tmpList = []
                tmpDict = {}
                inFile = open('miRNA/targetScan_miRNA_sets_geneSymbol.csv','r')
                inLines = [i.strip().split(',') for i in inFile.readlines() if i.strip()]
                for line in inLines:
                    if line[1].upper() in genesInBiclusters:
                        if not line[1] in tmpList:
                            tmpList.append(line[1].upper())
                        if not line[0] in tmpDict:
                            tmpDict[line[0]] = [line[1].upper()]
                        else:
                            tmpDict[line[0]].append(line[1].upper())
                inFile.close()
                pklFile = open('miRNA/targetScan.pkl','wb')
                cPickle.dump(tmpDict,pklFile)
                cPickle.dump(tmpList,pklFile)
            # Otherwise load the dumped pickle file if it exists
            else:
                print 'Loading pickled TargetScan predictions...'
                pklFile = open('miRNA/targetScan.pkl','rb')
                tmpDict = cPickle.load(pklFile)
                tmpList = cPickle.load(pklFile)
            pklFile.close()
            # Setup for analysis
            predDict = mgr.dict(tmpDict)
            pred_totalTargets = mgr.list()
            pred_totalTargets.append(set(tmpList))
            del tmpDict
            del tmpList
            print 'TargetScan has', len(predDict.keys()),'miRNAs.'

            def clusterHypergeo_ts(biclustId, db = predDict, allGenes = pred_totalTargets[0]):
                # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
                # Take gene list and compute overlap with each miRNA
                genes = allGenes.intersection(biclusters[biclustId].getGenes())
                writeMe = []
                keys1 = db.keys()
                m1s = []
                q = []
                m = []
                n = []
                k = []
                for m1 in keys1:
                    m1s.append(m1)
                    miRNAGenes = allGenes.intersection(db[m1])
                    q.append(len(set(miRNAGenes).intersection(genes)))
                    m.append(len(miRNAGenes))
                    n.append(len(allGenes)-len(miRNAGenes))
                    k.append(len(genes))
                results = phyper(q,m,n,k)
                min_miRNA = []
                perc_targets = []
                min_pValue = float(1)
                for i in range(1,len(results)):
                    if float(results[i]) <= float(0.05)/float(674) and not q[i]==0 and float(q[i])/float(k[i]) >= 0.1:
                        if min_miRNA==[] or float(results[i]) < min_pValue:
                            min_miRNA = [i]
                            perc_targets = [float(q[i])/float(k[i])]
                            min_pValue = float(results[i])
                        elif float(results[i])==min_pValue:
                            min_miRNA.append(i)
                            perc_targets.append(float(q[i])/float(k[i]))
                print 'Bicluster #', biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA])
                return [biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA]), ' '.join([str(j) for j in perc_targets]), min_pValue]

            # Set allGenes as intersect of pita_totalTargets and genesInBiclusters
            allGenes = pred_totalTargets[0].intersection(genesInBiclusters)

            # Run this using all cores available
            print 'Running TargetScan enrichment analyses...'
            cpus = cpu_count()
            biclustIds = biclusters.keys()
            pool = Pool(processes=cpus)
            res1 = pool.map(clusterHypergeo_ts,biclustIds)
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            pklFile = open('output/targetscan_3pUTR.pkl','wb')
            cPickle.dump(res1,pklFile)
        else:
            print 'Loading precached analysis...'
            pklFile = open('output/targetscan_3pUTR.pkl','rb')
            res1 = cPickle.load(pklFile)
        pklFile.close()

        # Stuff into biclusters
        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += miRNAInDict(m1.lower(),miRNAIDs)
            b1.addAttribute('targetscan_3pUTR',{'miRNA':r1[1], 'percentTargets':r1[2], 'pValue':r1[3], 'mature_seq_ids':miRNA_mature_seq_ids })
        print 'Done.\n'

        # TargetScan 3pUTR has been run on cMonkey run
        c1.targetscan_3pUTR = 1


    #################################################################
    ## Save out the final cMonkey object so we don't lose progress ##
    #################################################################
    if not os.path.exists('output/c1_all.pkl'):
        print 'Dumping Final cMonkey Object:'
        pklFile = open('output/c1_all.pkl','wb')
        cPickle.dump(c1,pklFile)
        pklFile.close()
        print 'Done.\n'


#################################################################
## Do postProcessing on cMonkey object                         ##
#################################################################
if not os.path.exists('output/c1_postProc.pkl'):
    #################################################################
    ## Save out the final cMonkey object so we don't lose progress ##
    #################################################################
    print 'Loading prechached cMonkey Object (c1_all.pkl):'
    pklFile = open('output/c1_all.pkl','rb')
    c1 = cPickle.load(pklFile)
    pklFile.close()
    print 'Done.\n'

    # Load up the expression ratios matrix
    ratioFile = gzip.open('../ratios.tsv.gz','rb')
    conditions = ratioFile.readline().split('\t')[1:]
    ratios = {}
    while 1:
        line = ratioFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        ratios[splitUp[0]] = dict(zip(conditions,splitUp[1:]))
    ratioFile.close()

    # Load the phenotype information
    # AGE,chemo_therapy,SURVIVAL,days_to_tumor_progression,SEX.bi,radiation_therapy,DEAD
    inFile = open('extras/phenotypes.csv','r')
    ids = inFile.readline().strip().split(',')[1:]
    phenotypes = {}
    for i in ids:
        phenotypes[i] = {}
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        phenotypes[splitUp[0]] = {}
        for i in range(len(ids)):
            phenotypes[ids[i]][splitUp[0]] = splitUp[i+1]
    inFile.close()

    def postProcess(bicluster):
        attributes = {}
        print ' Postprocessing cluster:', bicluster
        b1 = c1.getBicluster(bicluster)
        attributes['k'] = bicluster
        # Add number of genes and conditions
        attributes['k.rows'] = b1.getNumGenes()
        attributes['k.cols'] = b1.getNumConditions()
        # Get matrix of expression for genes
        genes = b1.getGenes()
        conditions = ratios[genes[0]].keys()
        matrix = [[ratios[gene][condition] for condition in conditions] for gene in genes]
        # Add first principal component variance explained
        pc1 = firstPrincipalComponent(matrix)
        attributes['pc1.var.exp'] = pc1[1]
        fpc = dict(zip(conditions, pc1[0]))
        attributes['pc1'] = fpc
        # Corrleation with patient traits
        cond2 = set(conditions).intersection(phenotypes['SURVIVAL'])
        pc1_1 = [fpc[i] for i in cond2]
        for phenotype in ['AGE','SEX.bi','chemo_therapy','radiation_therapy']:
            cond2 = set(conditions).intersection(phenotypes)
            p1_1 = [phenotypes[phenotype][i] for i in cond2]
            cor1 = correlation(pc1_1, p1_1)
            attributes[phenotype] = dict(zip(['rho','pValue'],cor1))
        # Association of bicluster expression with patient survival
        surv = [phenotypes['SURVIVAL'][i] for i in cond2]
        dead = [phenotypes['DEAD'][i] for i in cond2]
        age = [phenotypes['AGE'][i] for i in cond2]
        s1 = survival(surv, dead, pc1_1, age)
        attributes['Survival'] = dict(zip(['z','pValue'],s1[0]))
        attributes['Survival.AGE'] = dict(zip(['z','pValue'],s1[1]))
        return attributes

    # Do post processing
    print 'Do post processing...'
    cpus = cpu_count()
    biclustIds = c1.getBiclusters()
    pool = Pool(processes=cpus)
    res1 = pool.map(postProcess,biclustIds)
    pool.close()
    pool.join()
    print 'Done.\n'

    # Put results in cMonkey object
    for entry in res1:
        b1 = c1.getBicluster(entry['k'])
        for attribute in entry:
            if not attribute=='k':
                b1.addAttribute(attribute, entry[attribute])
    
    #################################################################
    ## TomTom Upstream motifs versus Jaspar and Transfac           ##
    #################################################################
    print 'Running TOMTOM on Upstream Motifs:'
    # Make needed directories
    if os.path.exists('tmp'):
        rmtree('tmp')
    if not os.path.exists('tmp/tomtom_out'):
        os.makedirs('tmp/tomtom_out')
    pssms = c1.getPssmsUpstream()
    upstreamMatches = {}
    if not os.path.exists('output/upstreamJasparTransfacComparison.pkl'):
        # Load JASPAR CORE Vertebarata motifs
        pklFile = open('motifs/jasparCoreVertebrata_redundant.pkl','rb')
        jasparPssms = cPickle.load(pklFile)
        for pssm1 in jasparPssms:
            jasparPssms[pssm1].setMethod('meme')
        pklFile.close()
        # Load Transfac 2012.1 motifs
        pklFile = open('motifs/transfac_2012.1_PSSMs_vertabrate.pkl','rb')
        transfacPssms = cPickle.load(pklFile)
        for pssm1 in transfacPssms:
            transfacPssms[pssm1].setMethod('meme')
        pklFile.close()
        # Uniprobe motifs
        pklFile = open('motifs/uniprobePSSMsNonRedundant.pkl','rb')
        uniprobePssms = cPickle.load(pklFile)
        for pssm1 in uniprobePssms:
            uniprobePssms[pssm1].setMethod('meme')
        pklFile.close()
        # Uniprobe motifs
        pklFile = open('motifs/selexPSSMsNonRedundant.pkl','rb')
        selexPssms = cPickle.load(pklFile)
        for pssm1 in selexPssms:
            selexPssms[pssm1].setMethod('meme')
        pklFile.close()

        # Write out results
        outFile = open('output/upstreamComparison_jaspar_transfac.csv','w')
        outFile.write('Motif Name,Original E-Value,Consensus,JASPAR Motif,JASPAR Consensus,TomTom.pValue,TomTom.qValue,Probe In Bicluster,Bicluster Residual')

        # Making MEME formatted files (makeFiles funciton in utils)
        print 'Making files...'
        makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=pssms.values(), targetPssms=jasparPssms.values(), num=1)
        makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=pssms.values(), targetPssms=transfacPssms.values(), num=2)
        makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=pssms.values(), targetPssms=uniprobePssms.values(), num=3)
        makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=pssms.values(), targetPssms=selexPssms.values(), num=4)

        # Run TomTom 
        print 'Comparing Upstream motifs against JASPAR CORE Vertebrata and Transfac 2012.1...'
        cpus = cpu_count()
        pool = Pool(processes=cpus)
        res1 = pool.map(runTomTom,[1,2,3,4])
        pool.close()
        pool.join()

        print 'Reading in Tomtom run...'
        outputFile = open('tmp/tomtom_out/tomtom1.out','r')
        outputFile.readline() # get rid of header
        outputLines = [ line.strip().split('\t') for line in outputFile.readlines() if float(line.split('\t')[4]) <= 0.05 ]
        outputFile.close()
        outputFile = open('tmp/tomtom_out/tomtom2.out','r')
        outputFile.readline() # get rid of header
        outputLines += [ line.strip().split('\t') for line in outputFile.readlines() if float(line.split('\t')[4]) <= 0.05 ]
        outputFile.close()
        outputFile = open('tmp/tomtom_out/tomtom3.out','r')
        outputFile.readline() # get rid of header
        outputLines += [ line.strip().split('\t') for line in outputFile.readlines() if float(line.split('\t')[4]) <= 0.05 ]
        outputFile.close()
        outputFile = open('tmp/tomtom_out/tomtom4.out','r')
        outputFile.readline() # get rid of header
        outputLines += [ line.strip().split('\t') for line in outputFile.readlines() if float(line.split('\t')[4]) <= 0.05 ]
        outputFile.close()

        # Now iterate through output and save data
        print 'Now parsing output for '+str(len(outputLines))+' matches...'
        for outputLine in outputLines:
            if len(outputLine)==9 and float(outputLine[4])<=0.05:
                tfName = outputLine[1].split('_')[0]
                if not outputLine[0] in upstreamMatches:
                    upstreamMatches[outputLine[0]] = [{'factor':outputLine[1],'confidence':outputLine[4]}]
                else:
                    upstreamMatches[outputLine[0]].append({'factor':outputLine[1],'confidence':outputLine[4]})
        outFile.close()
        pklFile = open('output/upstreamJasparTransfacComparison.pkl','wb')
        cPickle.dump(upstreamMatches,pklFile)
    else:
        print 'Loading precached upstream matches...'
        pklFile = open('output/upstreamJasparTransfacComparison.pkl','rb')
        upstreamMatches = cPickle.load(pklFile)
    pklFile.close()

    # Pump into pssms
    matched = []
    for pssmName in upstreamMatches:
        for match in upstreamMatches[pssmName]:
            pssms[pssmName].addMatch(factor=match['factor'], confidence=match['confidence'])
    print 'We matched '+str(len(upstreamMatches))+' upstream motifs.\n'


    #################################################################
    ## Get permuted p-values for upstream meme motifs              ##
    #################################################################
    # Make needed directories
    if os.path.exists('tmp'):
        rmtree('tmp')
    if not os.path.exists('tmp/tomtom_out'):
        os.makedirs('tmp/tomtom_out')
    # Compare the random motifs to the original motif in TOMTOM
    permPValues = {}
    matched = 0
    pssms = c1.getPssmsUpstream(de_novo_method='meme')
    if not os.path.exists('output/upstreamMotifPermutedPValues.csv'):
        outFile = open('output/upstreamMotifPermutedPValues.csv','w')
        outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
        pssmsNames = pssms.keys()
        print 'Loading precached random PSSMs...'
        randPssmsDict = {}
        for i in [5,10,15,20,25,30,35,40,45,50,55,60,65]:
            stdout.write(str(i)+' ')
            stdout.flush()
            pklFile = open(str(randPssmsDir)+'/pssms_upstream_'+str(i)+'.pkl','rb')
            randPssmsDict[i] = cPickle.load(pklFile)
            r1 = [randPssmsDict[i][pssm1].setMethod('meme') for pssm1 in randPssmsDict[i]]
            delMes = []
            for randPssm in randPssmsDict[i]:
                if not float(randPssmsDict[i][randPssm].getEValue()) <= float(maxEValue):
                    delMes.append(randPssm)
            for j in delMes:
                del randPssmsDict[i][j]

        print '\nMaking files...'
        for i in range(len(pssms)):
            clustSize = randPSSMClustSize((c1.getBicluster(int(pssmsNames[i].split('_')[0]))).getNumGenes())
            makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=[pssms[pssmsNames[i]]],targetPssms=randPssmsDict[clustSize].values(),num=i)

        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs avialable.'
        print 'Running TOMTOM to compare PSSMs...'
        pool = Pool(processes=cpus)
        pool.map(runTomTom,range(len(pssms)))
        pool.close()
        pool.join()

        print 'Reading in Tomtom run...'
        for run in range(len(pssms)):
            tomtomPValues = {}
            outputFile = open('tmp/tomtom_out/tomtom'+str(run)+'.out','r')
            output = outputFile.readlines()
            outputFile.close()
            # Now iterate through output and save data
            output.pop(0) # Get rid of header
            while len(output)>0:
                outputLine = output.pop(0).strip().split('\t')
                if len(outputLine)==9:
                    tomtomPValues[outputLine[1]] = float(outputLine[3])
            pValues = tomtomPValues.values()
            similar = 0
            for pValue in pValues:
                if float(pValue) <= float(0.05):
                    similar += 1
            # Write out the results
            mot = outputLine[0].split('_')[1]
            permPValues[outputLine[0]] = { mot+'.consensus':str(pssms[outputLine[0]].getConsensusMotif()), mot+'.permutedEV<=10':str(len(pValues)), mot+'.similar':str(similar), mot+'.permPV':str(float(similar)/float(1000)) }
            pssms[outputLine[0]].setPermutedPValue(str(float(similar)/float(1000)))
            outFile.write('\n'+str(outputLine[0])+',upstream,'+str(pssms[outputLine[0]].getEValue())+','+str(pssms[outputLine[0]].getConsensusMotif())+','+str(len(pValues))+','+str(similar)+','+str(1000)+','+str(float(similar)/float(1000)))
        outFile.close()
    else:
        print 'Using precalculated upstream p-values...'
        inFile = open('output/upstreamMotifPermutedPValues.csv','r')
        inFile.readline()
        upPValues = [i.strip().split(',') for i in inFile.readlines()]
        inFile.close()
        for line in upPValues:
            pssms[line[0]].setPermutedPValue(str(float(line[7])/float(1000)))
    print 'Done.\n'


    #################################################################
    ## Compare 3' UTR Weeder Motifs to miRBase using miRvestigator ##
    #################################################################
    print 'Running miRvestigator on 3\' UTR Motifs:'
    if not os.path.exists('output/m2m.pkl'):
        print 'Computing miRNA matches...'
        biclusters = c1.getBiclusters()
        pssms = c1.getPssms3pUTR()
        seqs3pUTR = c1.getSeqs3pUTR().values()
        m2m = miRvestigator(pssms.values(),seqs3pUTR,seedModel=[6,7,8],minor=True,p5=True,p3=True,wobble=False,wobbleCut=0.25)
        pklFile = open('output/m2m.pkl','wb')
        cPickle.dump(m2m,pklFile)
    else:
        print 'Loading precached miRNA matches...'
        pklFile = open('output/m2m.pkl','rb')
        m2m = cPickle.load(pklFile)
    pklFile.close()
    print 'Done.\n'


    #################################################################
    ## Convert miRNAs and Get permuted p-values for 3' UTR motifs  ##
    #################################################################
    pssms = c1.getPssms3pUTR()
    print 'Loading miRvestigator results...'
    # Convert miRvestigator results
    if not os.path.exists('output/miRvestigatorResults.pkl'):
        inFile = open('miRNA/scores.csv','r')
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]
        miRNA_matches = {}
        for line in lines:
            if not line[1]=='NA':
                miRNA_mature_seq_ids = []
                for i in line[1].split('_'):
                    miRNA_mature_seq_ids += miRNAInDict(i.lower(),miRNAIDs)
                miRNA_matches[line[0]] = {'miRNA':line[1],'model':line[2],'mature_seq_ids':miRNA_mature_seq_ids}
                for m1 in miRNA_mature_seq_ids:
                    pssms[line[0]].addMatch(factor=m1, confidence=line[2])
        pklFile = open('output/miRvestigatorResults.pkl','wb')
        cPickle.dump(miRNA_matches, pklFile)
    else:
        pklFile = open('output/miRvestigatorResults.pkl','rb')
        miRNA_matches = cPickle.load(pklFile)
        for m1 in miRNA_matches:
            for m2 in miRNA_matches[m1]['mature_seq_ids']:
                pssms[m1].addMatch(factor=m2, confidence=miRNA_matches[m1]['model'])
    pklFile.close()

    # Compile results to put them into the postProcessed
    print 'Get perumted p-values for 3\' UTR motifs...'
    pklFile = open('randPssms/weederRand.pkl','rb')
    weederRand = cPickle.load(pklFile)
    pklFile.close()
    pklFile = open('randPssms/weederRand_all.pkl','rb')
    weederRand_all = cPickle.load(pklFile)
    pklFile.close()
    clustSizes = sorted(weederRand['8bp'].keys())
    for pssm1 in pssms:
        seqNum = pssms[pssm1].getNumGenes()
        consensus = pssms[pssm1].getConsensusMotif()
        width = len(consensus)
        splitUp = pssm1.split('_')
        clustInd = 5
        for c in clustSizes:
            if seqNum > c:
                clustInd = c
        pValue = float(sum(1 for i in weederRand[str(width)+'bp'][clustInd] if float(i) >= float(pssms[pssm1].getEValue())))/float(len(weederRand[str(width)+'bp'][clustInd]))
        pValue_all = float(sum(1 for i in weederRand_all[str(width)+'bp'][clustInd] if float(i) >= float(pssms[pssm1].getEValue())))/float(len(weederRand_all[str(width)+'bp'][clustInd]))
        pssms[pssm1].setPermutedPValue({'pValue':pValue,'pValue_all':pValue_all})
    print 'Done.\n'


    #################################################################
    ## Run replication p-values                                    ##
    #################################################################
    # Dump a file containing all the genes for each cluster
    cmgFile = open('output/cluster.members.genes.txt','w')
    writeMe = []
    for b1 in c1.getBiclusters():
        writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).getGenes()))
    cmgFile.write('\n'.join(writeMe))
    cmgFile.close()

    # Dump a file containing all the genes for each cluster
    cmcFile = open('output/cluster.members.conditions.txt','w')
    writeMe = []
    for b1 in c1.getBiclusters():
        writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).getConditions()))
    cmcFile.write('\n'.join(writeMe))
    cmcFile.close()


    def runReplicaiton(repScript):
        # Fire up and run the replication
        repProc = Popen('cd '+repScript[0]+'; R --no-save < ' + repScript[1], shell=True, stdout=PIPE, stderr=PIPE)
        output = repProc.communicate()[0].split('\n')

    # Run replication on French dataset
    if (not os.path.exists('replication_French/replicationPvalues.csv')) or (not os.path.exists('replication_REMBRANDT/replicationPvalues.csv')) or (not os.path.exists('replication_GSE7696/replicationPvalues.csv')):
        print 'Run here..'
        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs avialable.'
        print 'Running TOMTOM to compare PSSMs...'
        pool = Pool(processes=cpus)
        pool.map(runReplication,[['replication_French','replicationDatasetPermutation.R'],['replication_REMBRANDT','replicationDatasetPermutation.R'],['replication_GSE7696','replicationDatasetPermutation.R']])
        pool.close()
        pool.join()

    #################################################################
    ## Read in replication p-values                                ##
    #################################################################
    # Read in replication p-values - French Dataset      
    # '','n.rows','overlap.rows','new.resid.norm.gbm','avg.norm.perm.resid.gbm','norm.perm.p.gbm','new.resid.norm.all','avg.norm.perm.resid.all','norm.perm.p.all','pc1.var.exp.gbm','avg.pc1.var.exp.gbm','pc1.perm.p.gbm','pc1.var.exp.all','avg.pc1.var.exp.all','pc1.perm.p.all','survival.gbm','survival.p.gbm','survival.age.gbm','survival.age.p.gbm','survival.all','survival.p.all','survival.age.all','survival.age.p.all'
    print 'Loading repliation p-values...'
    inFile = open('replication_French/replicationPvalues.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].replace('"','')))
        b1.addAttribute(key='replication_French',value={'French_new.resid.norm':splitUp[3], 'French_avg.resid.norm':splitUp[4], 'French_norm.perm.p':splitUp[5], 'French_pc1.var.exp':splitUp[9], 'French_avg.pc1.var.exp':splitUp[10], 'French_pc1.perm.p':splitUp[11], 'French_survival':splitUp[15], 'French_survival.p':splitUp[16], 'French_survival.age':splitUp[17], 'French_survival.age.p':splitUp[18]})
        b1.addAttribute(key='replication_French_all',value={'French_all_new.resid.norm':splitUp[6], 'French_all_avg.resid.norm':splitUp[7], 'French_all_norm.perm.p':splitUp[8], 'French_all_pc1.var.exp':splitUp[12], 'French_all_avg.pc1.var.exp':splitUp[13], 'French_all_pc1.perm.p':splitUp[14], 'French_all_survival':splitUp[19], 'French_all_survival.p':splitUp[20], 'French_all_survival.age':splitUp[21], 'French_all_survival.age.p':splitUp[22]})
    inFile.close()
    # Read in replication p-values - REMBRANDT Dataset      
    # "","n.rows","orig.resid","orig.resid.norm","overlap.rows","new.resid","avg.perm.resid","perm.p","new.resid.norm","avg.norm.perm.resid","norm.perm.p","survival","survival.p","survival.age","survival.age.p"
    print 'Loading repliation p-values...'
    inFile = open('replication_REMBRANDT/replicationPvalues.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].replace('"','')))
        b1.addAttribute(key='replication_REMBRANDT',value={'REMBRANDT_new.resid.norm':splitUp[3], 'REMBRANDT_avg.resid.norm':splitUp[4], 'REMBRANDT_norm.perm.p':splitUp[5], 'REMBRANDT_pc1.var.exp':splitUp[9], 'REMBRANDT_avg.pc1.var.exp':splitUp[10], 'REMBRANDT_pc1.perm.p':splitUp[11], 'REMBRANDT_survival':splitUp[15], 'REMBRANDT_survival.p':splitUp[16], 'REMBRANDT_survival.age':splitUp[17], 'REMBRANDT_survival.age.p':splitUp[18]})
        b1.addAttribute(key='replication_REMBRANDT_all',value={'REMBRANDT_all_new.resid.norm':splitUp[6], 'REMBRANDT_all_avg.resid.norm':splitUp[7], 'REMBRANDT_all_norm.perm.p':splitUp[8], 'REMBRANDT_all_pc1.var.exp':splitUp[12], 'REMBRANDT_all_avg.pc1.var.exp':splitUp[13], 'REMBRANDT_all_pc1.perm.p':splitUp[14], 'REMBRANDT_all_survival':splitUp[19], 'REMBRANDT_all_survival.p':splitUp[20], 'REMBRANDT_all_survival.age':splitUp[21], 'REMBRANDT_all_survival.age.p':splitUp[22]})
    # Read in replication p-values - GSE7696 Dataset
    # '', 'n.rows','overlap.rows','new.resid.norm.gbm','avg.norm.perm.resid.gbm','norm.perm.p.gbm','pc1.var.exp.gbm','avg.pc1.var.exp.gbm','pc1.perm.p.gbm','survival.gbm','survival.p.gbm','survival.age.gbm','survival.age.p.gbm'
    inFile = open('replication_GSE7696/replicationPvalues.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].replace('"','')))
        b1.addAttribute(key='replication_GSE7696',value={'GSE7696_new.resid.norm':splitUp[3], 'GSE7696_avg.resid.norm':splitUp[4], 'GSE7696_norm.perm.p':splitUp[5], 'GSE7696_pc1.var.exp':splitUp[6], 'GSE7696_avg.pc1.var.exp':splitUp[7], 'GSE7696_pc1.perm.p':splitUp[8], 'GSE7696_survival':splitUp[9], 'GSE7696_survival.p':splitUp[10], 'GSE7696_survival.age':splitUp[11], 'GSE7696_survival.age.p':splitUp[12]})
    inFile.close()

    ###########################################################################
    ## Run permuted p-value for variance epxlained first principal component ##
    ###########################################################################
    if not os.path.exists('output/residualPermutedPvalues_permAll.csv'):
        print 'Calculating FPC permuted p-values...'
        enrichProc = Popen("R --no-save < permutedResidualPvalues_permAll_mc.R", shell=True, stdout=PIPE, stderr=PIPE)
        output = enrichProc.communicate()[0]
        print 'Done.\n'

    #################################################################
    ## Read in residual permutations to use for filtering          ##
    #################################################################
    print 'Load residual permuted p-values...'
    inFile = open('output/residualPermutedPvalues_permAll.csv','r')
    # "","bicluster","n.rows","n.cols","orig.resid","avg.perm.resid","perm.p","orig.resid.norm","avg.norm.perm.resid","norm.perm.p","pc1.var.exp","avg.pc1.var.exp","pc1.perm.p"
    inFile.readline()
    inLines = inFile.readlines()
    for line in inLines:
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].strip('"')))
        b1.addAttribute(key='resid.norm.perm.p',value=str(splitUp[9]))
        b1.addAttribute(key='pc1.perm.p',value=str(splitUp[12]))
    inFile.close()
    print 'Done.\n'


    #################################################################
    ## Run functional enrichment and GO term similarity            ##
    #################################################################
    if not os.path.exists('funcEnrichment/biclusterEnrichment_GOBP.csv'):
        print 'Run functional enrichment...'
        enrichProc = Popen("cd funcEnrichment; R --no-save < enrichment.R", shell=True, stdout=PIPE, stderr=PIPE)
        output = enrichProc.communicate()[0]
        print 'Done.\n'
    if not os.path.exists('funcEnrichment/jiangConrath_hallmarks.csv'):
        print 'Run semantic similarity...'
        enrichProc = Popen("cd funcEnrichment; R --no-save < goSimHallmarksOfCancer.R", shell=True, stdout=PIPE, stderr=PIPE)
        output = enrichProc.communicate()[0]
        print 'Done.\n'

    #################################################################
    ## Read in hallmarks of cacner                                 ##
    #################################################################
    print 'Load Jiang-Conrath semantic similarity to Hallmarks of Cancer...'
    inFile = open('funcEnrichment/jiangConrath_hallmarks.csv','r')
    hallmarks = [i for i in inFile.readline().split(',') if not i.strip('"')=='']
    inLines = inFile.readlines()
    #hallmarksOfCancer = [int(line.split(',')[0].strip('"')) for line in inLines if [i for i in line.strip().split(',')[1:] if not i=='NA' and float(i) >= 0.8]]
    lines = [line.strip().split(',') for line in inLines]
    inFile.close()
    for line in lines:
        b1 = c1.getBicluster(int(line[0].strip('"')))
        b1.addAttribute(key='hallmarksOfCancer',value=dict(zip(hallmarks,line[1:])))
    print 'Done.\n'


    #################################################################
    ## Save out the final cMonkey object so we don't lose progress ##
    #################################################################
    print 'Dumping Final cMonkey Object:'
    pklFile = open('output/c1_postProc.pkl','wb')
    cPickle.dump(c1,pklFile)
    pklFile.close()
    print 'Done.\n'
else:
    print 'Loading from precached cMonkey Object:'
    pklFile = open('output/c1_postProc.pkl','rb')
    c1 = cPickle.load(pklFile)
    pklFile.close()
    print 'Done.\n'

#################################################################
## Write out the final post-processed file                     ##
#################################################################
print 'Write postProcessedVFinal.csv...'
postOut = []
hallmarksOfCancer = c1.getBicluster(1).getAttribute('hallmarksOfCancer').keys()
for i in sorted(c1.getBiclusters().keys()):
    writeMe = []
    b1 = c1.getBicluster(i)
    # Write file line by line
    #   a. Bicluster basics:  id, genes, conditions, resid, resid.norm, resid.norm.perm.p
    writeMe += [str(i), # id
                str(b1.getAttribute('k.rows')), # genes
                str(b1.getAttribute('k.cols')), # conditions
                str(b1.getNormResidual()), # normalized residual
                str(b1.getAttribute('resid.norm.perm.p')), # normalized residual permuted p-value
                str(b1.getAttribute('pc1.var.exp')), # Varience explained by first principal component
                str(b1.getAttribute('pc1.perm.p'))] # Varience explained by first principal component
    #   b. Upstream motifs:  meme.motif1.E, meme.motif1.consensus, meme.motif1.matches, meme.motif1.permPV
    motifNames = b1.getPssmsNamesUpstream()
    upstreamMotifs = { 'meme_motif1':None, 'meme_motif2':None, 'weeder_motif1':None, 'weeder_motif2':None }
    for m1 in motifNames:
        splitUp = m1.split('_')
        if splitUp[1]=='motif1' and splitUp[2]=='meme':
            upstreamMotifs['meme_motif1'] = m1
        if splitUp[1]=='motif2' and splitUp[2]=='meme':
            upstreamMotifs['meme_motif2'] = m1
        if splitUp[1]=='motif1' and splitUp[2]=='weeder':
            upstreamMotifs['weeder_motif1'] = m1
        if splitUp[1]=='motif2' and splitUp[2]=='weeder':
            upstreamMotifs['weeder_motif2'] = m1
    #   - MEME motifs
    for meme1 in ['meme_motif1','meme_motif2']:
        if not upstreamMotifs[meme1]==None:
            pssm1 = b1.getPssmUpstream(upstreamMotifs[meme1])
            matches = 'NA'
            if not pssm1.getMatches()==None:
                matches = ' '.join([match1['factor'] for match1 in pssm1.getMatches() if float(match1['confidence'])<=float(0.05)])
            writeMe += [str(pssm1.getEValue()), # E-value
                        str(pssm1.getPermutedPValue()), # Permuted p-value for motif
                        pssm1.getConsensusMotif(), # Motif consensus sequence
                        matches] # Matches to the motif from TransFac and Jaspar
        else:
            writeMe += ['NA', # E-value
                        'NA', # Permuted p-value for motif
                        'NA', # Motif consensus sequence
                        'NA'] # Matches to the motif from TransFac and Jaspar
    #   - WEEDER motifs
    for weeder1 in ['weeder_motif1','weeder_motif2']:
        if not upstreamMotifs[weeder1]==None:
            pssm1 = b1.getPssmUpstream(upstreamMotifs[weeder1])
            matches = 'NA'
            if not pssm1.getMatches()==None:
                matches = ' '.join([match1['factor'] for match1 in pssm1.getMatches() if float(match1['confidence'])<=float(0.05)])
            writeMe += [str(pssm1.getEValue()), # E-value
                        #str(pssm1.getPermutedPValue()), # Permuted p-value for motif
                        pssm1.getConsensusMotif(), # Motif consensus sequence
                        matches] # Matches to the motif from TransFac and Jaspar
        else:
            writeMe += ['NA', # E-value
                        #'NA', # Permuted p-value for motif
                        'NA', # Motif consensus sequence
                        'NA'] # Matches to the motif from TransFac and Jaspar

    #   c. 3' UTR motifs:  weeder.motif1.E, weeder.motif1.permPV, weeder.motif1.permPV_all, weeder.motif1.consensus, weeder.motif1.matches, weeder.motif1.model
    motifNames = b1.getPssmsNames3pUTR()
    p3utrMotifs = { 'weeder_motif1':None, 'weeder_motif2':None }
    for m1 in motifNames:
        splitUp = m1.split('_')
        if splitUp[1]=='motif1' and splitUp[2]=='weeder':
            p3utrMotifs['weeder_motif1'] = m1
        if splitUp[1]=='motif2' and splitUp[2]=='weeder':
            p3utrMotifs['weeder_motif2'] = m1
    #   - WEEDER motifs
    for weeder1 in ['weeder_motif1','weeder_motif2']:
        if not p3utrMotifs[weeder1]==None:
            pssm1 = b1.getPssm3pUTR(p3utrMotifs[weeder1])
            matches = 'NA'
            model = 'NA'
            if not pssm1.getMatches()==None:
                matches = ' '.join([miRNAIDs_rev[match1['factor']] for match1 in pssm1.getMatches()])
                model = pssm1.getMatches()[0]['confidence']
            permutedPValue = 'NA'
            if not pssm1.getPermutedPValue()==None:
                permutedPValue = pssm1.getPermutedPValue()
            writeMe += [str(pssm1.getEValue()), # E-value
                        str(permutedPValue['pValue']), # Permuted p-value for motif
                        str(permutedPValue['pValue_all']), # Permuted p-value for motif (all)
                        pssm1.getConsensusMotif(), # Motif consensus sequence
                        matches, # Matches to the motif to miRBase
                        model] # Model fit by miRvestigator
        else:
            writeMe += ['NA', # E-value
                        'NA', # Permuted p-value for motif
                        'NA', # Permuted p-value for motif (all)
                        'NA', # Motif consensus sequence
                        'NA', # Matches to the motif to miRBase
                        'NA'] # Model fit by miRvestigator

    #   d. Enriched miRNAs:  3pUTR_pita.miRNAs,3pUTR_pita.percTargets,3pUTR_pita.pValue,3pUTR_targetScan.miRNAs,3pUTR_targetScan.percTargets,3pUTR_targetScan.pValue
    for association in ['pita_3pUTR', 'targetscan_3pUTR']:
        a1 = b1.getAttribute(association)
        if not a1['miRNA']=='':
            writeMe += [str(a1['miRNA']).replace(';',' '), str(a1['percentTargets']).replace(';',' '), str(a1['pValue'])]
        else:
            writeMe += ['NA','NA','NA']

    #   e. Associations with traits:  age, sex.bi, chemo_therapy, radiation_therapy
    for association in ['AGE','SEX.bi', 'chemo_therapy','radiation_therapy']:
        ass1 = b1.getAttribute(association)
        writeMe += [str(ass1['rho']), str(ass1['pValue'])]
    surv1 = b1.getAttribute('Survival') 
    survAge1 = b1.getAttribute('Survival.AGE')
    writeMe += [str(surv1['z']), str(surv1['pValue']), str(survAge1['z']), str(survAge1['pValue'])]

    #   f. Replications:  'REMBRANDT_new.resid.norm','REMBRANDT_avg.resid.norm','REMBRANDT_norm.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p'
    replications_French = b1.getAttribute('replication_French')
    replications_REMBRANDT = b1.getAttribute('replication_REMBRANDT')
    replications_French_all = b1.getAttribute('replication_French_all')
    replications_REMBRANDT_all = b1.getAttribute('replication_REMBRANDT_all')
    replications_GSE7696 = b1.getAttribute('replication_GSE7696')
    for replication in ['French_pc1.var.exp','French_avg.pc1.var.exp','French_pc1.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p']:
        writeMe.append(str(replications_French[replication]))
    for replication in ['French_all_pc1.var.exp','French_all_avg.pc1.var.exp','French_all_pc1.perm.p','French_all_survival','French_all_survival.p','French_all_survival.age','French_all_survival.age.p']:
        writeMe.append(str(replications_French_all[replication]))
    for replication in ['REMBRANDT_pc1.var.exp','REMBRANDT_avg.pc1.var.exp','REMBRANDT_pc1.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p']:
        writeMe.append(str(replications_REMBRANDT[replication]))
    for replication in ['REMBRANDT_all_pc1.var.exp','REMBRANDT_all_avg.pc1.var.exp','REMBRANDT_all_pc1.perm.p','REMBRANDT_all_survival','REMBRANDT_all_survival.p','REMBRANDT_all_survival.age','REMBRANDT_all_survival.age.p']:
        writeMe.append(str(replications_REMBRANDT_all[replication]))
    for replication in ['GSE7696_pc1.var.exp','GSE7696_avg.pc1.var.exp','GSE7696_pc1.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p']:
        writeMe.append(str(replications_GSE7696[replication]))

    #   g. Glioma sub-group overlaps:  'NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM'
    #for overlap in ['NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']:
    #    writeMe.append(str(b1.getAttribute(overlap)))

    #   h. Hallmarks of Cancer:  Hanahan and Weinberg, 2011
    bhc1 = b1.getAttribute('hallmarksOfCancer')
    for hallmark in hallmarksOfCancer:
        writeMe.append(str(bhc1[hallmark]))

    # Add to the final output file
    postOut.append(deepcopy(writeMe))

postFinal = open('postProcessedVFinal.csv','w')
header = ['Bicluster', 'Genes', 'Patients', 'Norm. Residual', 'Norm. Residual Perm. P-Value', 'Var. Exp. First PC', 'Var. Exp. First PC Perm. P-Value', 'MEME Motif1 E-Value', 'MEME Motif1 Perm. P-Value', 'Up.MEME Motif1 Consensus', 'Up.MEME Motif1 Matches', 'Up.MEME Motif2 E-Value', 'Up.MEME Motif2 Perm. P-Value', 'Up.MEME Motif2 Consensus', 'Up.MEME Motif2 Matches', 'Up.WEEDER Motif1 Score', 'Up.WEEDER Motif1 Consensus', 'Up.WEEDER Motif1 Matches', 'Up.WEEDER Motif2 Score', 'Up.WEEDER Motif2 Consensus', 'Up.WEEDER Motif2 Matches', '3pUTR.WEEDER Motif1 E-Value', '3pUTR.WEEDER Motif1 Perm. P-Value', '3pUTR.WEEDER Motif1 Perm. P-Value (All)', '3pUTR.WEEDER Motif1 Consensus', '3pUTR.WEEDER Motif1 Matches', '3pUTR.WEEDER Motif1 Model', '3pUTR.WEEDER Motif2 E-Value', '3pUTR.WEEDER Motif2 Perm. P-Value', '3pUTR.WEEDER Motif2 Perm. P-Value (All)', '3pUTR.WEEDER Motif2 Consensus', '3pUTR.WEEDER Motif2 Matches', '3pUTR.WEEDER Motif2 Model', '3pUTR_pita.miRNAs', '3pUTR_pita.percTargets', '3pUTR_pita.pValue', '3pUTR_targetScan.miRNAs', '3pUTR_targetScan.percTargets', '3pUTR_targetScan.pValue', 'Age', 'Age.p', 'Sex', 'Sex.p', 'Chemotherapy', 'Chemotherapy.p', 'RadiationTherapy', 'RadiationTherapy.p', 'Survival', 'Survival.p', 'Survial.covAge', 'Survival.covAge.p', 'French_pc1.var.exp','French_avg.pc1.var.exp','French_pc1.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p', 'French_all_pc1.var.exp','French_all_avg.pc1.var.exp','French_all_pc1.perm.p','French_all_survival','French_all_survival.p','French_all_survival.age','French_all_survival.age.p', 'REMBRANDT_pc1.var.exp','REMBRANDT_avg.pc1.var.exp','REMBRANDT_pc1.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p', 'REMBRANDT_all_pc1.var.exp','REMBRANDT_all_avg.pc1.var.exp','REMBRANDT_all_pc1.perm.p','REMBRANDT_all_survival','REMBRANDT_all_survival.p','REMBRANDT_all_survival.age','REMBRANDT_all_survival.age.p', 'GSE7696_pc1.var.exp','GSE7696_avg.pc1.var.exp','GSE7696_pc1.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p'] + [i.strip() for i in hallmarksOfCancer]
postFinal.write(','.join(header)+'\n'+'\n'.join([','.join(i) for i in postOut]))
postFinal.close()
print 'Done.\n'

"""
names = ['id', 'k.rows', 'k.cols', 'resid', 'resid.norm', 'resid.norm.perm.p', 'motif1.E', 'motif1.consensus', 'motif1.matches', 'motif1.permutedEV<=10', 'motif1.permPV', 'motif2.E', 'motif2.consensus', 'motif2.matches', 'motif2.permutedEV<=10', 'motif2.permPV', '3pUTRmotif1.weederScore', '3pUTRmotif1.localPermP', '3pUTRmotif1.allPermP', '3pUTRmotif1.consensus', '3pUTRmotif1.miRNAs', '3pUTRmotif1.model', '3pUTR_pita.miRNAs', '3pUTR_pita.percTargets', '3pUTR_pita.pValue', '3pUTR_targetScan.miRNAs', '3pUTR_targetScan.percTargets', '3pUTR_targetScan.pValue', 'SEX.bi', 'SEX.bi.p', 'AGE', 'AGE.p', 'Survival', 'Survival.p', 'Survival.AGE', 'Survival.AGE.p', 'Survival.var', 'Survival.var.p', 'Survival.var.AGE', 'Survival.var.AGE.p','REMBRANDT_new.resid.norm','REMBRANDT_avg.resid.norm','REMBRANDT_norm.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p','NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']
names2 = ['k.rows', 'k.cols', 'resid', 'resid.norm', 'norm.perm.p', 'motif1.E', 'motif1.consensus', 'motif1.matches', 'motif1.permutedEV<=10', 'motif1.permPV', 'motif2.E', 'motif2.consensus', 'motif2.matches', 'motif2.permutedEV<=10', 'motif2.permPV', '3pUTRmotif1.eValue', '3pUTRmotif1.permPV.local', '3pUTRmotif1.permPV.all', '3pUTRmotif1.consensus', '3pUTRmotif1.miRNAs', '3pUTRmotif1.model', '3pUTR_pita.miRNAs', '3pUTR_pita.percTargets', '3pUTR_pita.pValue', '3pUTR_targetScan.miRNAs', '3pUTR_targetScan.percTargets', '3pUTR_targetScan.pValue', 'SEX.bi', 'SEX.bi.p', 'AGE', 'AGE.p', 'Survival', 'Survival.p', 'Survival.AGE', 'Survival.AGE.p', 'Survival.var', 'Survival.var.p', 'Survival.var.AGE', 'Survival.var.AGE.p','REMBRANDT_new.resid.norm','REMBRANDT_avg.resid.norm','REMBRANDT_norm.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p','NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']

#################################################################
## Making SIF File for Cytoscape                               ##
#################################################################
# Upstream motif comparisons
includedBiclusters = {}
pssmsUp = {}
pssms3p = {}
pita = {}
targetScan = {}
# Here is where I will filter based upon:
# 1. Survival <= ((0.05/720) = 6.94E-5)
# 2. Bicluster residual permutations <= ((0.05/720) = 6.94E-5)
# 3. A motif that fits one of the following:
#   a. Upstream motif with E-value <= 10 and Permuted P-value <= 0.05
#   b. 3' UTR motif with a perfect 8mer or 7mer match to miRBase miRNA seed sequence
inside = 0
m1 = 0
m2 = 0
wm = 0
p = 0
t = 0
for i in range(720):
    if (float(postProcessed[i+1]['Survival.AGE.p']) <= float(0.05)/float(720)) and (float(residPerms[i]) <= float(0.05)/float(720)) and (float(postProcessed[i+1]['resid.norm']) <= float(0.459)) and ((float(postProcessed[i+1]['REMBRANDT_norm.perm.p']) <= float(0.05)) or (float(postProcessed[i+1]['GSE7696_norm.perm.p']) <= float(0.05))): # and (i+1 in hallmarksOfCancer)
        inside += 1
        motif = 0
        if not postProcessed[i+1]['motif1.E']=='NA' and (((float(postProcessed[i+1]['motif1.E']) <= float(10)) and (float(postProcessed[i+1]['motif1.permPV']) <= 0.05)) or ((float(postProcessed[i+1]['motif1.E']) <= float(10)) and (not postProcessed[i+1]['motif1.matches']=='NA')) or ((float(postProcessed[i+1]['motif1.permPV']) <= float(0.05)) and (not postProcessed[i+1]['motif1.matches']=='NA'))):
            motif = 1
            m1 += 1
            pssmsUp[str(i+1)+'_motif1'] = c1.getBicluster(i+1).getPssmUpstream(str(i+1)+'_motif1')
        if not postProcessed[i+1]['motif2.E']=='NA' and (((float(postProcessed[i+1]['motif2.E']) <= float(10)) and (float(postProcessed[i+1]['motif2.permPV']) <= 0.05)) or ((float(postProcessed[i+1]['motif2.E']) <= float(10)) and (not postProcessed[i+1]['motif2.matches']=='NA')) or ((float(postProcessed[i+1]['motif2.permPV']) <= float(0.05)) and (not postProcessed[i+1]['motif2.matches']=='NA'))):
            motif = 1
            m2 += 1
            pssmsUp[str(i+1)+'_motif2'] = c1.getBicluster(i+1).getPssmUpstream(str(i+1)+'_motif2')
        if not postProcessed[i+1]['3pUTRmotif1.model']=='NA':
            motif = 1
            wm += 1
            pssms3p[str(i+1)+'_motif1'] = c1.getBicluster(i+1).getPssm3pUTR(str(i+1)+'_motif1')
        if ((not postProcessed[i+1]['3pUTR_pita.miRNAs']=='NA') and (float(postProcessed[i+1]['3pUTR_pita.percTargets'])>=float(0.5))):
            motif = 1
            p += 1
            pita[str(i+1)+'_pita'] = postProcessed[i+1]['3pUTR_pita.miRNAs']
        if ((not postProcessed[i+1]['3pUTR_targetScan.miRNAs']=='NA') and (float(postProcessed[i+1]['3pUTR_targetScan.percTargets'])>=float(0.5))):
            motif = 1
            t += 1
            targetScan[str(i+1)+'_targetScan'] = postProcessed[i+1]['3pUTR_targetScan.miRNAs']
        if motif==1:
            includedBiclusters[i+1] = c1.getBicluster(i+1)

print inside, m1, m2, wm, p, t, len(includedBiclusters)


# Identify similarity between upstream motifs using TomTom
pValueThreshold = 0.05/((len(pssmsUp)**2)/2)
tomtomUp = tomtom(pssmsUp.values(),pssmsUp.values(),c1.getNucFreqsUpstream(),'+ -',minOverlap=6)
putMeUp = dict(zip(pssmsUp.keys(),range(len(pssmsUp))))
pValues = []
pssmPssmPairs = []
for i in range(len(pssmsUp)):
    pValues.append(range(len(pssmsUp)))
for i in pssmsUp:
    for j in pssmsUp:
        pValue = tomtomUp.getScore(i,j)['pValue']

        pValues[putMeUp[i]][putMeUp[j]] = pValue
        if not i==j and float(pValue) <= float(pValueThreshold):
            pssmPssmPairs.append([i+'_Up','mm',j+'_Up'])
        if not i==j and not [i+'_Up','mm',j+'_Up'] in pssmPssmPairs and not postProcessed[int(i.split('_')[0])][i.split('_')[1]+'.matches']=='NA' and not postProcessed[int(j.split('_')[0])][j.split('_')[1]+'.matches']=='NA':
            matched1 = 0
            m1 = postProcessed[int(i.split('_')[0])][i.split('_')[1]+'.matches'].split(' ')
            m2 = postProcessed[int(j.split('_')[0])][j.split('_')[1]+'.matches'].split(' ')
            for match1 in m1:
                for match2 in m2:
                    if match1==match2:
                        pssmPssmPairs.append([i+'_Up','mm',j+'_Up'])
                        matched1 = 1
                        break
                if matched1 == 1:
                    break

# Also wouldn't hurt to add some connections between those motifs we identify via JASPAR and Transfac


# Now write out the matrix to a csv file for import into R
outFile = open('tmp/tomtom_pValuesUpstream.csv','w')
pssmsNames = pssmsUp.keys()
outFile.write(','+','.join(pssmsUp.keys()))
for i in range(len(pValues)):
    outFile.write('\n'+str(pssmsNames[i])+','+','.join(pValues[i]))
outFile.close()
if not os.path.exists('imgs'):
    os.mkdir('imgs')

# Plot PSSMs for all motifs
pssmsUp = c1.getPssmsUpstream()
for i in pssmsUp:
    pssmsUp[i].plot('imgs/all_'+i.split('_')[0]+'_'+pssmsUp[i].getConsensusMotif()+'.png')

pssmsUTR = c1.getPssms3pUTR()
for i in pssmsUTR:
    pssmsUTR[i].plot('imgs/all_'+i.split('_')[0]+'_'+pssmsUTR[i].getConsensusMotif()+'.png')


# 3' UTR motif comparisons. Instead of using TOMTOM use miRvestigator
p3UTRpssms = []
links = {}
motifNames = []
for i in all_miRNA_matches.keys():
    if i in pssms3p:
        motifNames.append(i+'_3pUTR')
    elif i in pita:
        motifNames.append(i)
    elif i in targetScan:
        motifNames.append(i)
for i in range(len(motifNames)):
    for j in range(i+1,len(motifNames)):
            doit = 0
            if motifNames[i].split('_')[1]=='targetScan':
                doit += 1
            if motifNames[i].split('_')[1]=='pita':
                doit += 1
            if motifNames[j].split('_')[1]=='targetScan':
                doit += 1
            if motifNames[j].split('_')[1]=='pita':
                doit += 1
            if doit<=1 and len(all_miRNA_matches[motifNames[i].rstrip('_3pUTR')].intersection(all_miRNA_matches[motifNames[j].rstrip('_3pUTR')])) > 0:
                m1 = ''
                m2 = ''
                if motifNames[i].split('_')[1]=='pita':
                    m1 = pita[motifNames[i]]
                elif motifNames[i].split('_')[1]=='targetScan':
                    m1 = targetScan[motifNames[i]]
                else:
                    m1 = motifNames[i]
                if motifNames[j].split('_')[1]=='pita':
                    m1 = pita[motifNames[j]]
                elif motifNames[j].split('_')[1]=='targetScan':
                    m1 = targetScan[motifNames[j]]
                else:
                    m2 = motifNames[j]
                pssmPssmPairs.append([m1,'mm',m2])

# First make gene to bicluster links
geneBiclusterPairs = []
biclusterPssmPairs = []
biclusterMiRNAPairs = []
biclusterClinicalTraits = []
biclusterClinicalTraitCorrelations = []
biclusterClinicalTraitPvalues = []
allGenes = []
negBiclusters = []
posBiclusters = []
allPssmsUp = []
allPssms3pUTR = []
allPita = []
allTargetScan = []
clinicalTraits = []
for i in includedBiclusters.keys():
    bi = c1.getBicluster(i)
    #if float(bi.getNormResidual())<=float(maxResidual):
    #if float(bi.getScore())<=float(maxScore) and float(bi.getSurvival()['"Survival"']['pValue'])<=float(maxSurv) and float(bi.getNormResidual())<=float(maxResidual):
    if float(bi.getSurvival()['"Survival.Age"']['zScore'])<0:
        negBiclusters.append(bi.getName().replace(' ','_'))
    else:
        posBiclusters.append(bi.getName().replace(' ','_'))
    # allBiclusters.append(bi.getName().replace(' ','_'))
    # Then make bicluster to pssm links
    tmpUp = bi.getPssmsUpstream()
    for pssm in tmpUp:
        if pssm.getName() in pssmsUp:
            allPssmsUp.append(pssm.getName()+'_Up')
            biclusterPssmPairs.append([bi.getName().replace(' ','_'),'bm',pssm.getName()+'_Up'])
    tmp3p = bi.getPssms3pUTR()
    for pssm in tmp3p:
        if pssm.getName() in pssms3p:
            allPssms3pUTR.append(pssm.getName()+'_3pUTR')
            biclusterPssmPairs.append([bi.getName().replace(' ','_'),'bm',pssm.getName()+'_3pUTR'])
    mirna_regulator = bi.getName().split(' ')[1]+'_pita'
    if bi.getName().split(' ')[1]+'_pita' in pita.keys():
        allPita.append(pita[mirna_regulator])
        biclusterPssmPairs.append([bi.getName().replace(' ','_'),'bm',pita[mirna_regulator]])
    mirna_regulator = bi.getName().split(' ')[1]+'_targetScan'
    if bi.getName().split(' ')[1]+'_targetScan' in targetScan.keys():
        allTargetScan.append(targetScan[mirna_regulator])
        biclusterPssmPairs.append([bi.getName().replace(' ','_'),'bm',targetScan[mirna_regulator]])
    genes = bi.getGenes()
    for gene in genes:
        allGenes.append(gene)
        geneBiclusterPairs.append([gene,'gb',bi.getName().replace(' ','_')])
    # Then add the clinical traits
    cor1 = bi.getCorrelations()
    clinicalTraits = cor1.keys()
    for trait in cor1:
        biclusterClinicalTraits.append([bi.getName().replace(' ','_'),'bt',trait])
        biclusterClinicalTraitCorrelations.append([bi.getName().replace(' ','_'),cor1[trait]['cor'],trait])
        if float(cor1[trait]['pValue'])==float(0):
            pValue = 50 # This means the p-value is 0, which can't happen. But it is below R's capabilities to compute.
        else:
            pValue = -log10(float(cor1[trait]['pValue']))
        biclusterClinicalTraitPvalues.append([bi.getName().replace(' ','_'),str(pValue),trait])
    
# Write SIF file
sifFile = open('cMonkey.sif','w')
sifFile.write('\n'.join(['\t'.join(nEn) for nEn in geneBiclusterPairs]))
sifFile.write('\n'+'\n'.join(['\t'.join(nEn) for nEn in biclusterPssmPairs]))
sifFile.write('\n'+'\n'.join(['\t'.join(nEn) for nEn in uniquify(pssmPssmPairs)]))
#sifFile.write('\n'+'\n'.join(['\t'.join(nEn) for nEn in uniquify(biclusterClinicalTraits)]))
sifFile.close()

# Write node attribute file
nodeAttFile = open('cMonkey_nodeAtt.txt','w')
nodeAttFile.write('geneBiclustPssm (class=Double)')
# nodeAttFile.write('Node Type')
allGenes = [[gene,' = 1'] for gene in allGenes]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allGenes]))
posBiclusters = [[bi,' = 2'] for bi in posBiclusters]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in posBiclusters]))
negBiclusters = [[bi,' = 3'] for bi in negBiclusters]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in negBiclusters]))
allPssmsUp1 = [[pssm,' = 4'] for pssm in allPssmsUp]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPssmsUp1]))
allPssms3pUTR1 = [[pssm,' = 5'] for pssm in allPssms3pUTR]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPssms3pUTR1]))
allPita = [[pssm,' = 6'] for pssm in allPita]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPita]))
allTargetScan = [[pssm,' = 7'] for pssm in allTargetScan]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allTargetScan]))
#clinicalTraits = [[trait,' = 5'] for trait in clinicalTraits]
#nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in clinicalTraits]))
nodeAttFile.close()

# Make node attribute file to add imgs for motifs
nodeAttFile = open('cMonkey_nodeAtt_imgs.txt','w')
nodeAttFile.write('imgs (class=string)')
allPssmsUp1 = [[pssm,' = file:///C:/Users/cplaisie/Desktop/glioma_final/network/final/imgs/',(str(pssm)+'.png')] for pssm in allPssmsUp]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPssmsUp1]))
allPssms3pUTR1 = [[pssm,' = file:///C:/Users/cplaisie/Desktop/glioma_final/network/final/imgs/',(str(pssm)+'.png')] for pssm in allPssms3pUTR]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPssms3pUTR1]))
nodeAttFile.close()

# Write node attribute file
edgeAttFile = open('cMonkey_edgeAtt_corClinicalTraits.txt','w')
edgeAttFile.write('corelationClinicalTraits (class=Double)')
writeMe = [str(edge[0])+' (bt) '+str(edge[2])+' = '+str(edge[1]) for edge in biclusterClinicalTraitCorrelations]
edgeAttFile.write('\n'+'\n'.join(writeMe))
edgeAttFile.close()

edgeAttFile = open('cMonkey_edgeAtt_pValClinicalTraits.txt','w')
edgeAttFile.write('pValueClinicalTraits (class=Double)')
writeMe = [str(edge[0])+' (bt) '+str(edge[2])+' = '+str(edge[1]) for edge in biclusterClinicalTraitPvalues]
edgeAttFile.write('\n'+'\n'.join(writeMe))
edgeAttFile.close()
"""
