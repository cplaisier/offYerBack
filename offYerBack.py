#################################################################
# @Program: offYerBack.py                                       #
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

from math import log10
import cPickle, os
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
from shutil import rmtree

### Parameters for filtering
maxScore = 0
maxResidual = 0.5 # 10
maxEValue = 10 # 1000000000000
maxSurv = 0.05/720
#pValueThreshold = 0.05/((375**2)/2)
postProcessedFile = 'postProcessed_REM.csv'
randPssmsDir = '/home/cplaisie/motifSimulations/10786_700bp/randPSSMs'

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

# Returns a list of all p-values
def getPValues(tomtomResults):
    pValues = []
    allScores = tomtomResults.getAllScores()
    o1 = allScores.keys()[0]
    for i in allScores[o1]:
        pValues.append(allScores[o1][i]['pValue'])
    return pValues

def runTomtom(i):
    runTomTom(i, distMeth='ed', qThresh='1', minOverlap=6) #blic5

# Make the files for a TomTom run
def makeFiles(nucFreqs, queryPssms, targetPssms, num, strands='+ -'):
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 3.0\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: '+strands+'\n\n'
    memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
    memeHeader += 'A '+str(round(float(nucFreqs['A']),3))+' C '+str(round(float(nucFreqs['C']),3))+' G '+str(round(float(nucFreqs['G']),3))+' T '+str(round(float(nucFreqs['T']),3))
    #extras = ['R','Y','K','M','S','W','B','D','H','V','N']
    #for e in extras:
    #    memeHeader += ' '+str(e)+' 0.000'
    # Make query PSSM file
    queryFile = open('tmp/query'+str(num)+'.tomtom','w')
    queryFile.write(memeHeader)
    for pssm in queryPssms:
        queryFile.write('\n\n'+pssm.getMemeFormatted())
    queryFile.close()
    # Make target PSSM file
    targetFile = open('tmp/target'+str(num)+'.tomtom','w')
    targetFile.write(memeHeader)
    for pssm in targetPssms:
        targetFile.write('\n\n'+pssm.getMemeFormatted())
    targetFile.close()

# Run TomTom on the files
def runTomTom(num, distMeth='ed', qThresh='1', minOverlap=6):
    # Arguments for tomtom
    tomtomArgs = ' -query tmp/query'+str(num)+'.tomtom -target tmp/target'+str(num)+'.tomtom -dist '+str(distMeth)+' -o tmp/tomtom_out -text -q-thresh '+str(qThresh)+' -min-overlap '+str(minOverlap)+' -verbosity 0'
    print tomtomArgs
    #p = Popen("tomtom" + tomtomArgs, shell=True)
    #sts = os.waitpid(p.pid, 0)
    errOut = open('tmp/stderr.out','w')
    tomtomProc = Popen("tomtom" + tomtomArgs, shell=True,stdout=PIPE, stderr=errOut)
    outputFile = open('tmp/tomtom_out/tomtom'+str(num)+'.out','w')
    output = tomtomProc.communicate()[0]
    outputFile.write(output)
    outputFile.close()
    errOut.close()

# Works!!!
from cMonkeyWrapper import cMonkeyWrapper
if not os.path.exists('c1.pkl'):
    c1 = cMonkeyWrapper('iter3000.RData')
    pklFile = open('c1.pkl','wb')
    cPickle.dump(c1,pklFile)
else:
    print 'Loading precached cMonkey object..'
    pklFile = open('c1.pkl','rb')
    c1 = cPickle.load(pklFile)
    print 'Done.'
pklFile.close()

"""
# Data from biclusters
print c1.getBiclusterNames()
print c1.getBicluster(1)
print biclusters[1]
print biclusters[1].getName()
print biclusters[1].getNumGenes()
print biclusters[1].getGenes()
print biclusters[1].getNumConditions()
print biclusters[1].getConditions()
print biclusters[1].getResidual()
print biclusters[1].getNormResidual()
print biclusters[1].getCorrelations()
print biclusters[1].getSurvival()
pssmsUp = biclusters[1].getPssmsUpstream()
print pssmsUp[0]
print pssmsUp[0].getMatrix()
print pssmsUp[0].getConsensusMotif()
print pssmsUp[0].getMemeFormatted()
seqUp = c1.getSeqsUpstream()
print seqUp[seqUp.keys()[0]]
pssmsUTR = biclusters[1].getPssms3pUTR()
print pssmsUTR[0]
print pssmsUTR[0].getMatrix()
print pssmsUTR[0].getConsensusMotif()
print pssmsUTR[0].getMemeFormatted()
seqUTR = c1.getSeqs3pUTR()
print seqUTR[seqUTR.keys()[0]]
"""
#####################
# Do postProcessing #
#####################
print 'Do post processing...'
stderrFile = open('stderr.out','w')
rProcess = Popen('R --no-save',shell=True,stdin=PIPE,stdout=PIPE,stderr=stderrFile)
stdout_val = rProcess.communicate('source(\'postProcessing.R\')\n')[0]
stderrFile.close()
print 'Done.'

# Read in results
# id, k.rows, k.cols, resid, resid.norm, motif1.E, motif2.E, 3pUtrMotif1.E, 3pUtrMotif2.E, SEX.bi, SEX.bi, AGE, AGE, GRADE.NUM, GRADE.NUM, MODEL.1, MODEL.1, SURVIVAL, SURVIVAL, ASTRO, ASTRO, GBM, GBM, CANCER, CANCER, Survival, Survival.p, Survival.AGE, Survival.AGE.p
postProcessed = {}
inFile = open('postProcessed_REM.csv','r')
postProcessedHeader = inFile.readline().strip()
names = ['id', 'k.rows', 'k.cols', 'resid', 'resid.norm', 'motif1.E', 'motif2.E', '3pUTRmotif1.E', '3pUTRmotif2.E', 'SEX.bi', 'SEX.bi.p', 'AGE', 'AGE.p', 'GRADE.NUM', 'GRADE.NUM.p', 'MODEL.1', 'MODEL.1.p', 'SURVIVAL', 'SURVIVAL.p', 'ASTRO', 'ASTRO.p', 'GBM', 'GBM.p', 'CANCER', 'CANCER.p', 'Survival', 'Survival.p', 'Survival.AGE', 'Survival.AGE.p']
for line in inFile.readlines():
    splitUp = line.strip().split(',')
    postProcessed[int(splitUp[0].lstrip('"').rstrip('"'))] = dict(zip(names,splitUp))
inFile.close()

###################################
## Compare miRNA motifs to miRDB ##
###################################
biclusters = c1.getBiclusters()
from miRNA2motifHMM import miRNA2motifHMM
pssms = (c1.getPssms3pUTR(maxSurv=maxSurv, maxNormResid=maxResidual)).values()
print 'Working on',len(pssms),'clusters:'
#pssms = [biclusters[313].getPssm3pUTR('313_motif1'), biclusters[55].getPssm3pUTR('55_motif2')]
seqs3pUTR = c1.getSeqs3pUTR().values()
m2m = miRNA2motifHMM(pssms,seqs3pUTR,0,8)
#print m2m.getTopHit('313_motif1')
#print m2m.getTopHit('55_motif2')

######################################
## Get permuted p-values for motifs ##
######################################
# Make needed directories
rmtree('tmp')
if not os.path.exists('tmp'):
    os.mkdir('tmp')
if not os.path.exists('tmp/tomtom_out'):
    os.mkdir('tmp/tomtom_out')

pssms = c1.getPssmsUpstream()
# Get permuted p-values for pssms
# Compare the random motifs to the original motif in TOMTOM
outFile = open('upstreamMotifPermutedPValues.csv','w')
outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
pssmsNames = pssms.keys()
print 'Making files...'
for i in range(len(pssms)):
    if not os.path.exists('tmp/query'+str(i)+'.tomtom') and not os.path.exists('tmp/target'+str(i)+'.tomtom'):
        clustSize = randPSSMClustSize((c1.getBicluster(int(pssmsNames[i].split('_')[0]))).getNumGenes())
        pklFile = open(str(randPssmsDir)+'/pssms_upstream_'+str(clustSize)+'.pkl','rb')
        randPssms = cPickle.load(pklFile)
        delMes = []
        for randPssm in randPssms:
            if not float(randPssms[randPssm].getEValue()) <= float(maxEValue):
                delMes.append(randPssm)
        for j in delMes:
            del randPssms[j]
        similar = 0
        if len(pssms)>0:
            makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=[pssms[pssmsNames[i]]],targetPssms=randPssms.values(),num=i)
print 'Done.'

# Run this using all cores available
cpus = cpu_count()
print 'There are', cpus,'CPUs avialable.' 
pool = Pool(processes=cpus)
pool.map(runTomtom,range(len(pssms)))
print 'Done with Tomtom runs.\n'

print 'Reading in Tomtom run...'
permPValues = {}
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
    outFile.write('\n'+str(outputLine[0])+',upstream,'+str(pssms[outputLine[0]].getEValue())+','+str(pssms[outputLine[0]].getConsensusMotif())+','+str(len(pValues))+','+str(similar)+','+str(1000)+','+str(float(similar)/float(1000)))
outFile.close()

# Now add them to the postProcessed file
for i in range(len(postProcessed)):
    if (str(i+1)+'_motif1') in permPValues:
        postProcessed[i+1] = dict(postProcessed[i+1], **permPValues[str(i+1)+'_motif1'])
    else:
        postProcessed[i+1] = dict(postProcessed[i+1], **{ 'motif1.consensus':'NA', 'motif1.permutedEV<=10':'NA', 'motif1.similar':'NA', 'motif1.permPV':'NA' })
    if (str(i+1)+'_motif2') in permPValues:
        postProcessed[i+1] = dict(postProcessed[i+1], **permPValues[str(i+1)+'_motif2'])
    else:
        postProcessed[i+1] = dict(postProcessed[i+1], **{ 'motif2.consensus':'NA', 'motif2.permutedEV<=10':'NA', 'motif2.similar':'NA', 'motif2.permPV':'NA' })

# Make needed directories
rmtree('tmp')
if not os.path.exists('tmp'):
    os.mkdir('tmp')
if not os.path.exists('tmp/tomtom_out'):
    os.mkdir('tmp/tomtom_out')

pssms = c1.getPssms3pUTR()
# Get permuted p-values for pssms
# Compare the random motifs to the original motif in TOMTOM
outFile = open('p3utrMotifPermutedPValues.csv','w')
outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
pssmsNames = pssms.keys()
print 'Making files...'
for i in range(len(pssms)):
    if not os.path.exists('tmp/query'+str(i)+'.tomtom') and not os.path.exists('tmp/target'+str(i)+'.tomtom'):
        clustSize = randPSSMClustSize((c1.getBicluster(int(pssmsNames[i].split('_')[0]))).getNumGenes())
        pklFile = open(str(randPssmsDir)+'/pssms_3pUTR_'+str(clustSize)+'.pkl','rb')
        randPssms = cPickle.load(pklFile)
        delMes = []
        for randPssm in randPssms:
            if not float(randPssms[randPssm].getEValue()) <= float(maxEValue):
                delMes.append(randPssm)
        for j in delMes:
            del randPssms[j]
        similar = 0
        if len(pssms)>0:
            makeFiles(nucFreqs=c1.getNucFreqs3pUTR(), queryPssms=[pssms[pssmsNames[i]]],targetPssms=randPssms.values(),num=i)
print 'Done.'

# Run this using all cores available
cpus = cpu_count()
print 'There are', cpus,'CPUs avialable.' 
pool = Pool(processes=cpus)
pool.map(runTomtom,range(len(pssms)))
print 'Done with Tomtom runs.\n'

print 'Reading in Tomtom run...'
permPValues = {}
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
    permPValues[outputLine[0]] = { '3pUTR'+mot+'.consensus':str(pssms[outputLine[0]].getConsensusMotif()), '3pUTR'+mot+'.permutedEV<=10':str(len(pValues)), '3pUTR'+mot+'.similar':str(similar), '3pUTR'+mot+'.permPV':str(float(similar)/float(1000)) }
    outFile.write('\n'+str(outputLine[0])+',upstream,'+str(pssms[outputLine[0]].getEValue())+','+str(pssms[outputLine[0]].getConsensusMotif())+','+str(len(pValues))+','+str(similar)+','+str(1000)+','+str(float(similar)/float(1000)))
outFile.close()

# Now add them to the postProcessed file
for i in range(len(postProcessed)):
    if (str(i+1)+'_motif1') in permPValues:
        postProcessed[i+1] = dict(postProcessed[i+1], **permPValues[str(i+1)+'_motif1'])
    else:
        postProcessed[i+1] = dict(postProcessed[i+1], **{ '3pUTRmotif1.consensus':'NA', '3pUTRmotif1.permutedEV<=10':'NA', '3pUTRmotif1.similar':'NA', '3pUTRmotif1.permPV':'NA' })
    if (str(i+1)+'_motif2') in permPValues:
        postProcessed[i+1] = dict(postProcessed[i+1], **permPValues[str(i+1)+'_motif2'])
    else:
        postProcessed[i+1] = dict(postProcessed[i+1], **{ '3pUTRmotif2.consensus':'NA', '3pUTRmotif2.permutedEV<=10':'NA', '3pUTRmotif2.similar':'NA', '3pUTRmotif2.permPV':'NA' })

# Then write out the final postprocessed file
names = ['id', 'k.rows', 'k.cols', 'resid', 'resid.norm', 'motif1.E', 'motif1.consensus', 'motif1.permutedEV<=10', 'motif1.permPV', 'motif2.E', 'motif2.consensus', 'motif2.permutedEV<=10', 'motif2.permPV', '3pUTRmotif1.E', '3pUTRmotif1.consensus', '3pUTRmotif1.permutedEV<=10', '3pUTRmotif1.permPV', '3pUTRmotif2.E', '3pUTRmotif2.consensus', '3pUTRmotif2.permutedEV<=10', '3pUTRmotif2.permPV', 'SEX.bi', 'SEX.bi.p', 'AGE', 'AGE.p', 'Survival', 'Survival.p', 'Survival.AGE', 'Survival.AGE.p']
names2 = ['k.rows', 'k.cols', 'resid', 'resid.norm', 'motif1.E', 'motif1.consensus', 'motif1.permutedEV<=10', 'motif1.permPV', 'motif2.E', 'motif2.consensus', 'motif2.permutedEV<=10', 'motif2.permPV', '3pUTRmotif1.E', '3pUTRmotif1.consensus', '3pUTRmotif1.permutedEV<=10', '3pUTRmotif1.permPV', '3pUTRmotif2.E', '3pUTRmotif2.consensus', '3pUTRmotif2.permutedEV<=10', '3pUTRmotif2.permPV', 'SEX.bi', 'SEX.bi.p', 'AGE', 'AGE.p', 'Survival', 'Survival.p', 'Survival.AGE', 'Survival.AGE.p']
postOut = []
for i in range(len(postProcessed)):
    postOut.append([str(i+1)])
    for nm1 in names2:
        postOut[i] += [postProcessed[i+1][nm1]]
postFinal = open('postProcessedV2.csv','w')
postFinal.write(','.join(names)+'\n'+'\n'.join([','.join(i) for i in postOut]))
postFinal.close()

#################################
# Making SIF File for Cytoscape #
#################################
pssms = c1.getPssmsUpstream(maxScore=maxScore,maxEValue=maxEValue,maxSurv=maxSurv, maxNormResid=maxResidual)
pValueThreshold = 0.05/((len(pssms)**2)/2)
from tomtom import tomtom
tomtomUp = tomtom(pssms.values(),pssms.values(),c1.getNucFreqsUpstream(),'+ -',minOverlap=6)
putMeUp = dict(zip(pssms.keys(),range(len(pssms))))
pValues = []
pssmPssmPairs = []
for i in range(len(pssms)):
    pValues.append(range(len(pssms)))
for i in pssms:
    for j in pssms:
        pValue = tomtomUp.getScore(i,j)['pValue']
        pValues[putMeUp[i]][putMeUp[j]] = pValue
        if not i==j and float(pValue) <= float(pValueThreshold):
            pssmPssmPairs.append([i+'_Up','mm',j+'_Up'])
# Now write out the matrix to a csv file for import into R
outFile = open('tmp/tomtom_pValuesUpstream.csv','w')
pssmsNames = pssms.keys()
outFile.write(','+','.join(pssms.keys()))
for i in range(len(pValues)):
    outFile.write('\n'+str(pssmsNames[i])+','+','.join(pValues[i]))
outFile.close()

# UTR Works!!!
pssms = c1.getPssms3pUTR(maxScore=maxScore,maxEValue=maxEValue,maxSurv=maxSurv, maxNormResid=maxResidual)
pValueThreshold = 0.05/((len(pssms)**2)/2)
tomtom3pUTR = tomtom(pssms.values(),pssms.values(),c1.getNucFreqs3pUTR(),'+')
putMe3pUTR = dict(zip(pssms.keys(),range(len(pssms))))
pValues = []
for i in range(len(pssms)):
    pValues.append(range(len(pssms)))
for i in pssms:
    for j in pssms:
        curScore = tomtom3pUTR.getScore(i,j)
        if curScore['orientation']=='+': #and float(curScore['pValue'])<=0.0000431034:
            pValues[putMe3pUTR[i]][putMe3pUTR[j]] = curScore['pValue']
            if not i==j and float(curScore['pValue'])<= float(pValueThreshold):
                pssmPssmPairs.append([i+'_3pUTR','mm',j+'_3pUTR'])
        else:
            pValues[putMe3pUTR[i]][putMe3pUTR[j]] = '1'
# Now write out the matrix to a csv file for import into R
outFile = open('tmp/tomtom_pValuesUTR.csv','w')
pssmsNames = pssms.keys()
outFile.write(','+','.join(pssms.keys()))
for i in range(len(pValues)):
    outFile.write('\n'+str(pssmsNames[i])+','+','.join(pValues[i]))
outFile.close()

# First make gene to bicluster links
geneBiclusterPairs = []
biclusterPssmPairs = []
biclusterClinicalTraits = []
biclusterClinicalTraitCorrelations = []
biclusterClinicalTraitPvalues = []
allGenes = []
negBiclusters = []
posBiclusters = []
allPssmsUp = []
allPssms3pUTR = []
clinicalTraits = []
for i in c1.getBiclusterNames():
    bi = c1.getBicluster(i)
    #if float(bi.getNormResidual())<=float(maxResidual):
    if float(bi.getScore())<=float(maxScore) and float(bi.getSurvival()['"Survival.Age"']['pValue'])<=float(maxSurv) and float(bi.getNormResidual())<=float(maxResidual):
        if float(bi.getSurvival()['"Survival.Age"']['zScore'])<0:
            negBiclusters.append(bi.getName().replace(' ','_'))
        else:
            posBiclusters.append(bi.getName().replace(' ','_'))
        # allBiclusters.append(bi.getName().replace(' ','_'))
        genes = bi.getGenes()
        for gene in genes:
            allGenes.append(gene)
            geneBiclusterPairs.append([gene,'gb',bi.getName().replace(' ','_')])
        # Then make bicluster to pssm links
        pssmsUp = bi.getPssmsUpstream()
        for pssm in pssmsUp:
            if float(pssm.getEValue())<=float(maxEValue):
                allPssmsUp.append(pssm.getName()+'_Up')
                biclusterPssmPairs.append([bi.getName().replace(' ','_'),'bm',pssm.getName()+'_Up'])
        pssms3pUTR = bi.getPssms3pUTR()
        for pssm in pssms3pUTR:
            if float(pssm.getEValue())<=float(maxEValue):
                allPssms3pUTR.append(pssm.getName()+'_3pUTR')
                biclusterPssmPairs.append([bi.getName().replace(' ','_'),'bm',pssm.getName()+'_3pUTR'])
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
allPssmsUp = [[pssm,' = 4'] for pssm in allPssmsUp]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPssmsUp]))
allPssms3pUTR = [[pssm,' = 5'] for pssm in allPssms3pUTR]
nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in allPssms3pUTR]))
#clinicalTraits = [[trait,' = 5'] for trait in clinicalTraits]
#nodeAttFile.write('\n'+'\n'.join([''.join(att) for att in clinicalTraits]))
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

