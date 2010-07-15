#################################################################
# @Program: permutationCoexpressionV4.py                        #
# @Version: 4                                                   #
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
# Copyrighted by Chris Plaisier  3/4/2010                       #
#################################################################

import os
from subprocess import *

# Parameters for permutation
clusterSizes = 'c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)'
pValueThreshold = 0.05
permutations = 1000

# Run mast
def mast(queryFile=None, seqFile=None, bgFile=None, ev=99999, mt=0.99):
    mastArgs = str(queryFile)+' '+str(seqFile)+' -bfile tmp/meme/bgFile.meme -nostatus -stdout -text -seqp -ev 9999999 -mev 9999999 -mt 1 -remcorr -brief'
    print mastArgs
    #errOut = open('tmp/meme/stderr.out','w')
    mastProc = Popen("mast " + mastArgs, shell=True,stdout=PIPE) #,stderr=errOut)
    output = mastProc.communicate()[0].split('\n')
    # Write out to a file
    outFile = open('out.txt','w')
    outFile.write('\n'.join(output))
    outFile.close()
    # Return results like cMonkey
    res1 = {}
    for i in range(len(output)):
        if 'COMBINED P-VALUE' in output[i]:
            splitUp = [j for j in output[i].strip().split(' ') if j]
            print splitUp
            res1[output[i-2].strip()] = { 'pValue': float(splitUp[6]), 'eValue': float(splitUp[9]) }
    return res1

# Run fimo
def fimo(queryFile=None, seqFile=None, outputPT=0.05):
    fimoArgs = '--bgfile tmp/meme/bgFile.meme --o tmp --output-pthresh '+str(outputPT)+' --text --verbosity 1 '+str(queryFile)+' '+str(seqFile)
    print fimoArgs
    #errOut = open('tmp/meme/stderr.out','w')
    fimoProc = Popen("~/bin/fimo " + fimoArgs, shell=True,stdout=PIPE) #,stderr=errOut)
    output = fimoProc.communicate()[0].split('\n')
    # Write out to a file
    outFile = open('out.txt','w')
    outFile.write('\n'.join(output))
    outFile.close()
    # Return results like cMonkey
    res1 = {}
    output.pop(0)
    # Motif	Seq	Start	Stop	Log-odds	p-value	Site
    for line1 in output:
        splitUp = line1.split('\t')
        if len(splitUp)==7:
            if not splitUp[1] in res1:
                res1[splitUp[1]] = { str(splitUp[2])+'_'+str(splitUp[3]):{ 'orientation':splitUp[0][0], 'start':splitUp[2], 'stop':splitUp[3], 'logOdds':float(splitUp[4]), 'pValue': float(splitUp[5]), 'site': splitUp[6] } }
            else:
                res1[splitUp[1]][str(splitUp[2])+'_'+str(splitUp[3])] = { 'orientation':splitUp[0][0], 'start':splitUp[2], 'stop':splitUp[3], 'logOdds':float(splitUp[4]), 'pValue': float(splitUp[5]), 'site': splitUp[6] }
    return res1

# Function to sort the bins by size
def qsortBasedOn(sortMe, basedOn):
    if not len(sortMe) == len(basedOn):
        return 'ERROR!'
    if len(basedOn) <= 1:
            return [sortMe, basedOn]
    pivot = basedOn.pop(0)
    pivotSM = sortMe.pop(0)
    greater = []
    lesser = []
    greaterSM = []
    lesserSM = []
    while len(basedOn) > 0:
        cur = basedOn.pop(0)
        curSM = sortMe.pop(0)
        if cur >= pivot:
            greater.append(cur)
            greaterSM.append(curSM)
        else:
            lesser.append(cur)
            lesserSM.append(curSM)
    greaterOut = qsortBasedOn(greaterSM, greater)
    lesserOut = qsortBasedOn(lesserSM, lesser)
    return [lesserOut[0] + [pivotSM] + greaterOut[0], lesserOut[1] + [pivot] + greaterOut[1]]

# Benjamini-Hochberg - takes a dictionary of { name: pValue, ... }
def benjaminiHochberg(dict1, tests=10768, alpha=0.05):
    # First sort the results
    sorted1 = qsortBasedOn(dict1.keys(), dict1.values())[0]
    # Then control based on FDR
    res1 = []
    alpha = float(alpha)
    #res1 = [sorted1[i] for i in range(len(sorted1)) if dict1[sorted1[i]] <= alpha/float(tests-i)]
    for i in range(len(sorted1)):
        if dict1[sorted1[i]] <= alpha/float(tests-i):
            res1.append(sorted1[i])
        else:
            break
    return res1

# Load up the cMonkey run to get permuted pValues
from cMonkeyWrapper import cMonkeyWrapper
c1 = cMonkeyWrapper('iter3000.RData')

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

# Build an R script to gather the data from the cMonkey run
sendToR = []
# First load up the cmonkey run
sendToR.append('load(\'iter3000.RData\')')
sendToR.append('library(cMonkey)')
# Get a subset
pssmsUp = { '469_motif1':pssmsUp['469_motif1'], '150_motif1':pssmsUp['150_motif1'], '237_motif1':pssmsUp['237_motif1'], '237_motif2':pssmsUp['237_motif2'], '427_motif1':pssmsUp['427_motif1'], '122_motif2':pssmsUp['122_motif2'], '616_motif1':pssmsUp['616_motif1'], '616_motif2':pssmsUp['616_motif2'], '401_motif1':pssmsUp['401_motif1'], '401_motif2':pssmsUp['401_motif2'], '333_motif1':pssmsUp['333_motif1'], '333_motif2':pssmsUp['333_motif2']}
#pssmsUp = { '237_motif1':pssmsUp['237_motif1'] }
seqs = c1.getSeqsUpstream()
# Make sequence database (fasta)
seqFile = open('tmp/seqDB.fa','w')
seqFile.write('\n'.join(['>'+seq+'\n'+seqs[seq] for seq in seqs]))
# Read in background frequencies
nucFreqs = c1.getNucFreqsUpstream()

# Load up all the pssms with an e-value less than 10
pssmsUp = c1.getPssmsUpstream(maxEValue=10)

# Load up all Biclusters
biclusters = c1.getBiclusters()
for bicluster in biclusters:
    # If motifs are decent then use them and make a page for them
    pssms = [i for i in [str(bicluster)+'_motif1', str(bicluster)+'_motif2'] if i in pssmsUp ]
    if len(pssms)>0:
        for pssm1 in pssms:
            # 1. Identify genes w/ motif out of all genes from a cMonkey run
            #   a. Read in motif (pssm object).
            #   b. Use mast to identify genes, and a cutoff to discretize.
            print 'In Silico predition of motif targets...'
            # Read in PSSMs
            #pssms = c1.getBicluster(427).getPssmsUpstream()
            print pssmsUp[pssm1].getConsensusMotif()
            pssms = [pssmsUp[pssm1]]
            
            ### Run FIMO ###
            fimoHeader = ''
            fimoHeader += 'MEME version 3.0\n\n'
            fimoHeader += 'ALPHABET= ACGT\n\n'
            # Here is where we tell it what strand: for miRNAs this would just be '+'
            fimoHeader += 'strands: + -\n\n'
            fimoHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
            fimoHeader += 'A '+str(round(float(nucFreqs['A']),3))+' C '+str(round(float(nucFreqs['C']),3))+' G '+str(round(float(nucFreqs['G']),3))+' T '+str(round(float(nucFreqs['T']),3))
            queryFile = open('tmp/query'+str(pssmsUp[pssm1].getName())+'.fimo','w')
            queryFile.write(fimoHeader)
            for pssm in pssms:
                queryFile.write('\n\n'+pssm.getMemeFormatted())
            queryFile.close()

            # Run FIMO - Benjamini-Hochberg corrected
            f1 = fimo(queryFile='tmp/query'+str(pssmsUp[pssm1].getName())+'.fimo', seqFile='tmp/seqDB.fa', outputPT=0.05)
            genes = []
            pValues = []
            for gene in f1:
                k1 = f1[gene].keys()
                genes.append(gene)
                if len(f1[gene])==1:
                    pValues.append(f1[gene][k1[0]]['pValue'])
                else:
                    min = f1[gene][k1.pop(0)]['pValue']
                    for inst1 in k1:
                        if float(f1[gene][inst1]['pValue']) < float(min):
                            min = f1[gene][inst1]['pValue']
                    pValues.append(min)
            
            # Determine which genes have motif
            haveMotif = benjaminiHochberg(dict(zip(genes,pValues)),tests=10768,alpha=0.05)
            outFile = open('motifSet'+str(pssmsUp[pssm1].getName())+'.txt','w')
            outFile.write('\n'.join(haveMotif))
            outFile.close()
            print 'In Silico Predicted Targets:',len(haveMotif),'out of',len(f1.keys())
            print 'Done.\n'
            
            # 2. Calculate residual for x number of permuations for each cluster size of genes w/ motif
            #   a. Write out gene set
            #   b. Fire up R
            #   c. Get a random sample
            #   d. Calculate residuals
            if(len(haveMotif) >= 70):
                print 'Doing permutation with genes w/ motifs...'
                # Read genes w/ motif
                sendToR.append('l1 = read.csv(\'motifSet'+str(pssmsUp[pssm1].getName())+'.txt\', header=F)')
                # Get random samples
                sendToR.append('clustRange = c(5, 10, 15, 20, 30, 40, 50, 60, 70)')
                sendToR.append('r1 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('r3 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('r2 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('r4 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('t1 = matrix(nrow=length(clustRange),ncol=1,dimnames=list(clustRange,\'T-Test P-Value\'))')
                sendToR.append('t2 = matrix(nrow=length(clustRange),ncol=1,dimnames=list(clustRange,\'T-Test P-Value\'))')
                sendToR.append('mu1 = c(1:length(clustRange))')
                sendToR.append('mu2 = c(1:length(clustRange))')
                sendToR.append('k.cols = e$get.cols('+str(pssmsUp[pssm1].getName().split('_')[0])+')')
                sendToR.append('for(i in 1:length(clustRange)) {')
                sendToR.append('    r1[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(l1[,1], clustRange[i]), k.cols, TRUE) } )')
                sendToR.append('    r2[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(rownames(e$ratios$ratios), clustRange[i]), k.cols, TRUE) } )')
                sendToR.append('    r3[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(l1[,1], clustRange[i]), colnames(e$ratios$ratios), TRUE) } )')
                sendToR.append('    r4[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(rownames(e$ratios$ratios), clustRange[i]), colnames(e$ratios$ratios), TRUE) } )')
                sendToR.append('    mu1[i] = mean(r1[,i])')
                sendToR.append('    mu2[i] = mean(r2[,i])')
                sendToR.append('    mu3[i] = mean(r3[,i])')
                sendToR.append('    mu4[i] = mean(r4[,i])')
                sendToR.append('    tmp1 = t.test(r1[,i],r2[,i])')
                sendToR.append('    t1[i] = tmp1$p.value*(tmp1$stat/abs(tmp1$stat))')
                sendToR.append('    tmp2 = t.test(r3[,i],r4[,i])')
                sendToR.append('    t2[i] = tmp2$p.value*(tmp2$stat/abs(tmp2$stat))')
                sendToR.append('}')
                sendToR.append('write.csv(r1,file=\'r1'+str(pssmsUp[pssm1].getName())+'.csv\')')
                sendToR.append('write.csv(r2,file=\'r2'+str(pssmsUp[pssm1].getName())+'.csv\')')
                sendToR.append('write.csv(r3,file=\'r3'+str(pssmsUp[pssm1].getName())+'.csv\')')
                sendToR.append('write.csv(r4,file=\'r4'+str(pssmsUp[pssm1].getName())+'.csv\')')
                sendToR.append('write.csv(data.frame(mu1=mu1, mu2=mu2, t1, mu3=mu3, mu4=mu4, t2),file=\'t/t_'+str(pssmsUp[pssm1].getName())+'.csv\')')
    # Now add in analysis for combined motifs if there are 2
    if len(pssms)==2:
        ### Run MAST ###
        # Make files - query file (list of log odds ratios)
        mastHeader = 'ALPHABET= ACGT\n'
        queryFile = open('tmp/query'+str(bicluster)+'.mast','w')
        queryFile.write(mastHeader)
        queryFile.write('\n'.join([pssm.getMastFormatted() for pssm in pssms]))
        queryFile.close()
        m1 = mast(queryFile='tmp/query'+str(bicluster)+'.mast', seqFile='tmp/seqDB.fa', bgFile='tmp/meme/bgFile.meme')
        
        # Determine which genes have motif
        haveMotif = []
        for g1 in m1.keys():
            if m1[g1]['pValue'] <= pValueThreshold:
                haveMotif.append(g1)
        outFile = open('motifSet'+str(bicluster)+'.txt','w')
        outFile.write('\n'.join(haveMotif))
        outFile.close()
        
        if(len(haveMotif) >= 70):
                print 'Doing permutation with genes w/ motifs...'
                # Read genes w/ motif
                sendToR.append('l1 = read.csv(\'motifSet'+str(bicluster)+'.txt\', header=F)')
                # Get random samples
                sendToR.append('clustRange = c(5, 10, 15, 20, 30, 40, 50, 60, 70)')
                sendToR.append('r1 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('r3 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('r2 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('r4 = matrix(ncol=length(clustRange),nrow='+str(permutations)+')')
                sendToR.append('t1 = matrix(nrow=length(clustRange),ncol=1,dimnames=list(clustRange,\'T-Test P-Value\'))')
                sendToR.append('t2 = matrix(nrow=length(clustRange),ncol=1,dimnames=list(clustRange,\'T-Test P-Value\'))')
                sendToR.append('mu1 = c(1:length(clustRange))')
                sendToR.append('mu2 = c(1:length(clustRange))')
                sendToR.append('k.rows = e$get.rows('+str(bicluster)+')')
                sendToR.append('k.cols = e$get.cols('+str(bicluster)+')')
                sendToR.append('o0 = e$residual.submatrix(e$ratios$ratios, k.rows, colnames(e$ratios$ratios))')
                sendToR.append('o1 = e$residual.submatrix(e$ratios$ratios, k.rows, colnames(e$ratios$ratios), TRUE)')
                sendToR.append('o2 = e$residual.submatrix(e$ratios$ratios, k.rows, k.cols)')
                sendToR.append('o3 = e$residual.submatrix(e$ratios$ratios, k.rows, k.cols, TRUE)')
                sendToR.append('o4 = e$cluster.resid('+str(bicluster)+')')
                sendToR.append('o5 = e$cluster.resid('+str(bicluster)+', varNorm=TRUE)')
                sendToR.append('for(i in 1:length(clustRange)) {')
                sendToR.append('    r1[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(l1[,1], clustRange[i]), k.cols, TRUE) } )')
                sendToR.append('    r2[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(rownames(e$ratios$ratios), clustRange[i]), k.cols, TRUE) } )')
                sendToR.append('    r3[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(l1[,1], clustRange[i]), colnames(e$ratios$ratios), TRUE) } )')
                sendToR.append('    r4[,i] = sapply(1:'+str(permutations)+', function(x) { e$residual.submatrix(e$ratios$ratios, sample(rownames(e$ratios$ratios), clustRange[i]), colnames(e$ratios$ratios), TRUE) } )')
                sendToR.append('    mu1[i] = mean(r1[,i])')
                sendToR.append('    mu2[i] = mean(r2[,i])')
                sendToR.append('    mu3[i] = mean(r3[,i])')
                sendToR.append('    mu4[i] = mean(r4[,i])')
                sendToR.append('    tmp1 = t.test(r1[,i],r2[,i])')
                sendToR.append('    t1[i] = tmp1$p.value*(tmp1$stat/abs(tmp1$stat))')
                sendToR.append('    tmp2 = t.test(r3[,i],r4[,i])')
                sendToR.append('    t2[i] = tmp2$p.value*(tmp2$stat/abs(tmp2$stat))')
                sendToR.append('}')
                sendToR.append('write.csv(r1,file=\'r1'+str(bicluster)+'.csv\')')
                sendToR.append('write.csv(r2,file=\'r2'+str(bicluster)+'.csv\')')
                sendToR.append('write.csv(r3,file=\'r3'+str(bicluster)+'.csv\')')
                sendToR.append('write.csv(r4,file=\'r4'+str(bicluster)+'.csv\')')
                sendToR.append('write.csv(data.frame(mu1=mu1, mu2=mu2, t1, mu3=mu3, mu4=mu4, t2),file=\'t/t_'+str(bicluster)+'.csv\')')
                sendToR.append('write.csv(rbind(c('all.conds','all.conds.norm','clust.conds','clust.conds.norm','clust.resid','clust.resid.norm'),c(o0,o1,o2,o3,o4,o5)),file=\'t/o'+str(bicluster)+'.csv\')')

######################
#!  Run the R code  !#
print 'Accessing RData object.'
stderrFile = open('stderr.out','w')
rProcess = Popen('R --no-save',shell=True,stdin=PIPE,stdout=PIPE,stderr=stderrFile)
stdout_val = rProcess.communicate('\n'.join(sendToR))[0]
stderrFile.close()
#!  Done!           !#
######################
print 'Done.\n'

