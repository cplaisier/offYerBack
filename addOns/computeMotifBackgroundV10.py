#################################################################
# @Program: computeMotifBackground.py                           #
# @Version: 10                                                  #
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
# Copyrighted by Chris Plaisier  2/11/2009                      #
#################################################################

from multiprocessing import Pool, cpu_count, Manager
from pssm import pssm
from copy import deepcopy
from subprocess import *
import cPickle, os
import matplotlib
matplotlib.use('Agg')
from pylab import hist, savefig, close, plot, title, xlabel, ylabel, xlim
from numpy import array, float64, log10
from weblogolib import *
import numpy, corebio
from time import gmtime, strftime

# Make needed directories
if not os.path.exists('tmp'):
    os.mkdir('tmp')
if not os.path.exists('tmp/tomtom_out'):
    os.mkdir('tmp/tomtom_out')

# Plot a PSSM using weblogo
def plotPssm(pssm, fileName):
    dist = numpy.array( pssm.getMatrix(), numpy.float64 ) 
    data = LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist*100)
    options = LogoOptions()
    options.color_scheme = colorscheme.nucleotide
    format = LogoFormat(data, options)
    fout = open(fileName, 'w')
    png_formatter(data, format, fout)
    fout.close()

# Plot a PSSM using weblogo
def plotPssmMP(inArr1):
    pssm = inArr1[0]
    fileName = inArr1[1]
    dist = numpy.array( pssm.getMatrix(), numpy.float64 ) 
    data = LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist*100)
    options = LogoOptions()
    options.color_scheme = colorscheme.nucleotide
    format = LogoFormat(data, options)
    fout = open(fileName, 'w')
    png_formatter(data, format, fout)
    fout.close()

def convertPSSM(pssmMatrix):
    output = ['\tPOSITION\tA\tC\tG\tT']
    for i in range(1,(len(pssmMatrix)+1)):
        output.append(str(i)+'\t'+str(pssmMatrix[i-1][0])+'\t'+str(pssmMatrix[i-1][1])+'\t'+str(pssmMatrix[i-1][2])+'\t'+str(pssmMatrix[i-1][3]))
    return output

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

# Get the maximum number of hits
def getMaxHit(pssmHits):
    max = 'None'
    for pssm in pssmHits:
        if pssmHits[pssm] > 0 and (max=='None' or pssmHits[pssm] > pssmHits[max]):
            max = pssm
    return max

def getPValue(pssm1, pssm2, tomtomPValues):
    if pssm1 in tomtomPValues:
        if pssm2 in tomtomPValues[pssm1]:
            return tomtomPValues[pssm1][pssm2]
    return tomtomPValues[pssm2][pssm1]

# A recursive functon to bin the motifs
def binMotifs(pssm1, pssms2, tomtomPValues, pssms, eValueThreshold, pValueThreshold, orig):
    bin = []
    pssms2.remove(pssm1)
    for pssm2 in pssms2:
        if float(pssms[pssm1].getEValue())<=float(eValueThreshold) and float(pssms[pssm2].getEValue())<=float(eValueThreshold) and float(getPValue(pssm1, pssm2, tomtomPValues)) <= float(pValueThreshold) and float(getPValue(orig, pssm2, tomtomPValues)) <= float(0.001):
            bin += binMotifs(pssm2, pssms2, tomtomPValues, pssms, eValueThreshold, pValueThreshold, orig)
    return [pssm1]+bin

# Function to sort the bins by size
def qsortBins(binned):
    if len(binned) <= 1:
            return binned
    pivot = binned.pop(0)
    greater = qsortBins([i for i in binned if len(i) > len(pivot)])
    lesser_eq = qsortBins([i for i in binned if len(i) <= len(pivot)])
    return greater + [pivot] + lesser_eq

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

# Function to sort the bins by size
def qsortBasedOnUp(sortMe, basedOn):
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
    return greaterOut[0] + [pivotSM] + lesserOut[0]

# Flips the nucleotide values for a pssm matrix [[A0, C0, G0, T0], ..., [An, Cn, Gn, Tn]]
def flip(pssmMatrix):
    flipped = []
    for i in pssmMatrix:
        flipped = [[i[3], i[2], i[1], i[0]]] + flipped
    return flipped

# Returns a merged pssm object for the mergedSet
def merge2PSSMs(name, pssm1, pssm2, weight1, weight2, offset, overlap, orientation):
    pssm1Mat = deepcopy(pssm1.getMatrix())
    # Flip PSSM matrix values for nucleotides if orientation=='-'
    pssm2Mat = deepcopy(pssm2.getMatrix())
    if orientation=='-':
        pssm2Mat = flip(pssm2Mat)
    # Merge PSSM matrices
    newPSSMMat = []
    newWeight = []
    if offset>0:
        for i in range(offset):
            row1 = pssm1Mat.pop(0)
            w1 = weight1.pop(0)
            newPSSMMat += [[row1[0]*100*w1, row1[1]*100*w1, row1[2]*100*w1, row1[3]*100*w1]]
            newWeight += [w1]
    elif offset<0:
        for i in range(-offset):
            row2 = pssm2Mat.pop(0)
            w2 = weight2.pop(0)
            newPSSMMat += [[row2[0]*100*w2, row2[1]*100*w2, row2[2]*100*w2, row2[3]*100*w2]]
            newWeight += [w2]
    maxLen = 0
    if len(pssm1Mat)>=len(pssm2Mat):
        maxLen = len(pssm1Mat)
    else:
        maxLen = len(pssm2Mat)
    for i in range(maxLen):
        if len(pssm1Mat)>0:
            row1 = pssm1Mat.pop(0)
            w1 = weight1.pop(0)
        else:
            row1 = [0, 0, 0, 0]
            w1 = 0
        if len(pssm2Mat)>0:
            row2 = pssm2Mat.pop(0)
            w2 = weight2.pop(0)
        else:
            row2 = [0, 0, 0, 0]
            w2 = 0
        newPSSMMat += [[row1[0]*100*w1+row2[0]*100*w2, row1[1]*100*w1+row2[1]*100*w2, row1[2]*100*w1+row2[2]*100*w2, row1[3]*100*w1+row2[3]*100*w2]]
        newWeight += [w1+w2]
    # Normalize matrix
    for i in range(len(newPSSMMat)):
        total = newPSSMMat[i][0]+newPSSMMat[i][1]+newPSSMMat[i][2]+newPSSMMat[i][3]
        newPSSMMat[i][0] = newPSSMMat[i][0]/total
        newPSSMMat[i][1] = newPSSMMat[i][1]/total
        newPSSMMat[i][2] = newPSSMMat[i][2]/total
        newPSSMMat[i][3] = newPSSMMat[i][3]/total
    # Now build the new PSSM object
    newPSSM = pssm(biclusterName=name, nsites='30', eValue='0', pssm=newPSSMMat, genes=[])
    print name, newPSSM.getConsensusMotif(), pssm1.getName(), pssm1.getConsensusMotif(), pssm2.getName(), pssm2.getConsensusMotif(), offset, overlap, orientation
    return [newPSSM, newWeight]

"""# Get complete linkage p-value
def completeLinkage(a, b, tomtomPValues, pValueThreshold):
    maxP = float(0)
    for i in a:
        for j in b:
            if tomtomPValues[i][j]['pValue']>maxP:
                maxP = tomtomPValues[i][j]['pValue']
                # If the maximum p-value exceeds our thrshold don't bother continuing
                if maxP>pValueThreshold:
                    return maxP
    return maxP
"""

# Get complete linkage p-value
def completeLinkage(a, b, tomtomPValues, pValueThreshold):
    d1 = []
    for i in a:
        for j in b:
            c = tomtomPValues[i][j]['pValue']
            if c > pValueThreshold:
                return c
            d1.append(c)
    return max(d1)

# Now read in nucleotide frequencies
nucFreqsReg = {}
nucFreqsFile = open('nucFreqs.csv','r')
nucFreqsFile.readline() # Skip header
upFreq = nucFreqsFile.readline().strip().split(',')
nucFreqsReg['upstream'] = {'A':upFreq[1],'C':upFreq[2],'G':upFreq[2],'T':upFreq[1]}
p3utrFreq = nucFreqsFile.readline().strip().split(',')
nucFreqsReg['3pUTR'] = {'A':p3utrFreq[1],'C':p3utrFreq[2],'G':p3utrFreq[2],'T':p3utrFreq[1]}
nucFreqsFile.close()

#regions = ['upstream','3pUTR']
regions = ['upstream']
clusterSizes = range(4, 65)
#clusterSizes = [30]
#clusterSizes = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
pValueThreshold = float(0.01)
eValueThreshold = float(10)
print 'P-Value Threshold = ',pValueThreshold, '; E-Value Threshold = ', eValueThreshold
# Build an html page for all the clusterSizes
if not os.path.exists('results'):
    os.mkdir('results')
overallHtml = open('results/index.html','w')
overallHtml.write('<html>\n<head><title>Summary of Randomized MEME Runs</title></head>\n<body style=\'font-family: arial, verdana, sans-serif\'>\n')
for region in regions:
    motifBins = {}
    # For each region plot the number of E-Values less than specified threshold
    goodEValsAll = []
    overallHtml.write('<center><h1>'+region+'</h1></center>\n<center><a href=\''+region+'/imgs/barPlot.png\'><img src=\''+region+'/imgs/barPlot.png\' height=200></a></center>\n')
    # Make a new table fore each region
    overallHtml.write('<center><table border=1 bgcolor=\'#336699\'><tr><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Cluster Size</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Histogram of E-Values</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Motifs (E-Values <= '+str(eValueThreshold)+')</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Motif Bins (P-Vaue <= '+str(pValueThreshold)+')</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Singletons</font></th></tr>\n')
    overallHtml.write('<tr><td colspan=5 bgcolor=\'#FFFFCC\'><center><b><a href=\''+region+'/superBin/clustering.html\'>SuperBin</a></b></center></td></tr>\n')
    if not os.path.exists('results/'+region):
        os.mkdir('results/'+region)
    for clusterSize in clusterSizes:
        motifBins[clusterSize] = {}
        print 'CLUSTER SIZE = '+str(clusterSize)
        # Start a shared memory manager
        mgr = Manager()
        nucFreqs = mgr.dict(nucFreqsReg[region])
        # Setup for the Tomtom runs
        pklFile = open('randPSSMs/pssms_'+str(region)+'_'+str(clusterSize)+'.pkl','rb')
        pssms = cPickle.load(pklFile)
        # Plot the histogram
        if not os.path.exists('results/'+region+'/imgs'):
            os.mkdir('results/'+region+'/imgs')
        eVals = []
        goodEVals = []
        for pssm1 in pssms:
            eVals.append(float(pssms[pssm1].getEValue()))
            if float(pssms[pssm1].getEValue()) <= eValueThreshold:
                goodEVals.append(pssm1)
        goodEValsAll.append(len(goodEVals))
        eVals = array(eVals)
        eVals = log10(eVals)
        h1 = hist(eVals,100)
        h1 = xlim(-20,10)
        savefig('results/'+region+'/imgs/hist'+str(clusterSize)+'.png',format='png')
        close()
        csDir = 'results/'+region+'/cs'+str(clusterSize)
        if not os.path.exists(csDir):
            os.mkdir(csDir)
        if not os.path.exists(csDir+'/imgs'):
            os.mkdir(csDir+'/imgs')
        overallHtml.write('<tr><td bgcolor=\'#FFFFCC\'><center><a href=\''+region+'/cs'+str(clusterSize)+'/clustering.html\'>'+str(clusterSize)+'</a></center></td><td bgcolor=\'#FFFFCC\'><center><a href=\''+region+'/imgs/hist'+str(clusterSize)+'.png\'><img src=\''+region+'/imgs/hist'+str(clusterSize)+'.png\' height=100></a></center></td><td bgcolor=\'#FFFFCC\'><center>'+str(len(goodEVals))+'</center></td>')
        # Back to loading stuff up
        if not os.path.exists('tomtomPValues_'+str(clusterSize)+'.pkl'):
            # First make all files where each query file is one pssm and the target file is all remaining pssms, do the upper triangle
            pssmNames = deepcopy(goodEVals)
            print 'Making Files ('+str(len(goodEVals))+')...'
            runs = 0
            for i in range(len(goodEVals)):
                queryPssms = [pssms[goodEVals[i]]]
                targetPssms = []
                for other in goodEVals:
                    targetPssms.append(pssms[other])
                makeFiles(nucFreqs=nucFreqs, queryPssms=queryPssms, targetPssms=targetPssms, num=i)
            print 'Done.'
            # Run this using all cores available
            cpus = cpu_count()
            print 'There are', cpus,'CPUs avialable.' 
            pool = Pool(processes=cpus)
            pool.map(runTomtom,range(len(goodEVals)))
            print 'Done with Tomtom runs.\n'
            tomtomPValues = {}
            print 'Reading in Tomtom run...'
            for run in range(len(goodEVals)):
                outputFile = open('tmp/tomtom_out/tomtom'+str(run)+'.out','r')
                output = outputFile.readlines()
                outputFile.close()
                # Now iterate through output and save data
                output.pop(0) # Get rid of header
                while len(output)>0:
                    outputLine = output.pop(0).strip().split('\t')
                    if len(outputLine)==9:
                        if not outputLine[0] in tomtomPValues:
                            tomtomPValues[outputLine[0]] = { outputLine[1]: { 'pValue': float(outputLine[3]), 'qValue':outputLine[4], 'offset':outputLine[2], 'overlap': outputLine[5], 'orientation': outputLine[8] } }
                        else:
                            tomtomPValues[outputLine[0]][outputLine[1]] = { 'pValue': float(outputLine[3]), 'qValue':outputLine[4], 'offset':outputLine[2], 'overlap': outputLine[5], 'orientation': outputLine[8] }
                        if outputLine[1] in tomtomPValues and outputLine[0] in tomtomPValues[outputLine[1]]:
                            if tomtomPValues[outputLine[0]][outputLine[1]]['pValue'] > tomtomPValues[outputLine[1]][outputLine[0]]['pValue']:
                                tomtomPValues[outputLine[0]][outputLine[1]] = tomtomPValues[outputLine[1]][outputLine[0]]
                            else:
                                tomtomPValues[outputLine[0]][outputLine[1]] = tomtomPValues[outputLine[1]][outputLine[0]]
            tomtomFile = open('tomtomPValues_'+str(clusterSize)+'.pkl','wb')
            cPickle.dump(tomtomPValues,tomtomFile)
        else:
            print 'Reading in Tomtom run...'
            tomtomFile = open('tomtomPValues_'+str(clusterSize)+'.pkl','rb')
            tomtomPValues = cPickle.load(tomtomFile)
            # Should be removed later
            for p1 in tomtomPValues:
                for p2 in tomtomPValues[p1]:
                    tomtomPValues[p1][p2]['pValue'] = float(tomtomPValues[p1][p2]['pValue'])
        tomtomFile.close()
        print 'Done.\n'

        ## Agglomerative complete linkage clustering algorithm
        #   Step 1. Find minimum pairwise PSSM p-value between all clusters
        #   Step 2. Combine clusters
        #   Step 3. Continue until p-value between clusters is greater than pValueThreshold
        print 'Agglomerative complete linkage clustering start...'
        clusters = [[i] for i in deepcopy(goodEVals)]
        print strftime("%a, %d %b %Y %H:%M:%S", gmtime())
        while 1:
            # Find minimum pairwise PSSM p-value (complete linkage)
            clusters2 = deepcopy(clusters)
            minP = ['', '', float(1)]
            for c1 in clusters:
                clusters2.pop(0)
                for c2 in clusters2:
                    testMe = completeLinkage(c1, c2, tomtomPValues, pValueThreshold)
                    #testMe = max([max(v1) for v1 in [[tomtomPValues[p1][p2]['pValue'] for p2 in c2] for p1 in c1]])
                    if testMe<minP[2]:
                        minP = [c1, c2, testMe]
            if minP[2]<pValueThreshold:
                clusters.remove(minP[0])
                clusters.remove(minP[1])
                clusters.append(minP[0]+minP[1])
            else:
                break
        print len(clusters), clusters
        print strftime("%a, %d %b %Y %H:%M:%S", gmtime())
        print 'Done.\n'

        # Remove singletons and sort bins
        print 'Sort bins and remove singletons...'
        binned = []
        binSimP = {}
        singletons = []
        for bin1 in clusters:
            if len(bin1)>1:
                # Pick the exemplar motif
                if len(bin1)>2:
                    # Pick the most representative PSSM
                    pssmHits = {}
                    for pssm1 in bin1:
                        pssmHits[pssm1] = 0
                        for pssm2 in bin1:
                            if not pssm1==pssm2 and tomtomPValues[pssm1][pssm2]['pValue'] <= pValueThreshold:
                                pssmHits[pssm1] += 1
                    maxHit = getMaxHit(pssmHits)
                    tmpBin = deepcopy(bin1)
                    tmpBin.remove(maxHit)
                    # Sort the remaining based on similarity to maximum hit
                    basedOn = []
                    for i in tmpBin:
                        basedOn.append(tomtomPValues[maxHit][i]['pValue'])
                    sortedOut = qsortBasedOn(tmpBin,basedOn)
                    binned.append(deepcopy([maxHit]+sortedOut[0]))
                    binSimP.update(dict(zip([maxHit]+sortedOut[0],['Reference']+sortedOut[1])))
                    motifBins[clusterSize][bin1[0]] = { 'pssm': pssms[maxHit], 'numHits': len(bin1) }
                else:
                    motifBins[clusterSize][bin1[0]] = { 'pssm': pssms[bin1[0]], 'numHits': len(bin1) }
                    binned.append(deepcopy(bin1))
                    binSimP.update(dict(zip(bin1,['-','-'])))
            else:   
                singletons.append(bin1[0])
        binned = qsortBins(binned)
        print 'Done.\n'
    
        overallHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+str(len(binned))+'</center></td><td bgcolor=\'#FFFFCC\'><center>'+str(len((singletons)))+'</center></td></tr>\n')

        # Write out the bins
        print 'Write out bins...'
        clusteringHtml = open(csDir+'/clustering.html','w')
        clusteringHtml.write('<html>\n<head><title>Clustering of MEME Results for Cluster Size '+str(clusterSize)+'</title>\n<style type=\'text/css\'>\n.hidden {\ndisplay: none;\n }\n</style>\n</head>\n<body style=\'font-family: arial, verdana, sans-serif\'>\n')
        clusteringHtml.write('<center><p><h2>Cluster Size = '+str(clusterSize)+'</h2></p></center><center><table border=1 bgcolor=\'#336699\'><tr><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Merged Name</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Merged PSSM</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Consensus</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Number PSSMs Merged</font></th></tr>\n')
        clusteringHtmlEnd = ''
        print 'Plot pssms...'
        cpus = cpu_count()
        pool = Pool(processes=cpus)
        # Build the array for plotting PSSMs [[pssm, fileName], ... ]
        inArr1 = []
        for pssm1 in goodEVals:
            inArr1.append([pssms[pssm1],csDir+'/imgs/'+str(pssm1)+'.png'])
        pool.map(plotPssmMP,inArr1)
        print 'Done plotting.'
        for bin1 in binned:
            #outFile.write('\n'+str(repPssm)+','+pssms[repPssm].getConsensusMotif()+','+str(len(binnedDict[repPssm]))+','+' '.join(bin))
            clusteringHtml.write('<tr>')
            # Plot the maxHit PSSM
            clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+bin1[0]+'</center></td>') # Best hit
            pssmImgName = 'imgs/'+str(bin1[0])+'.png'
            #plotPssm(pssms[bin1[0]], csDir+'/'+pssmImgName)
            clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><img src=\''+pssmImgName+'\'></td>') # Best hit PSSM Plot
            clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+pssms[bin1[0]].getConsensusMotif()+'</center></td>') # Best hit Consensus
            clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center><b><a href=\'bin'+bin1[0]+'.html\'>'+str(len(bin1))+'</a></b></center></td>') # Link to see all motifs in cluster
            clusteringHtml.write('</tr>\n')
            # Added so we can send the best hit PSSM to STAMP
            clusteringHtmlEnd += '<div id=\''+str(bin1[0])+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
            clusteringHtmlEnd += '<span class=\'gaggle-name hidden\'>'+str(bin1[0])+'</span>\n'
            clusteringHtmlEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
            clusteringHtmlEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[bin1[0]].getMatrix()))+'x4</span>\n'
            clusteringHtmlEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
            clusteringHtmlEnd += '\n'.join(convertPSSM(pssms[bin1[0]].getMatrix()))
            clusteringHtmlEnd += '</div>\n</div>\n'
            # Now build each clusters html file
            clusterHtml = open(csDir+'/bin'+bin1[0]+'.html','w')
            clusterHtml.write('<html>\n<head><title>Bin '+bin1[0]+' of MEME Results for Cluster Size '+str(clusterSize)+'</title>\n<style type=\'text/css\'>\n.hidden {\ndisplay: none;\n }\n</style>\n</head>\n<body style=\'font-family: arial, verdana, sans-serif\'>\n')
            clusterHtml.write('<center><table border=1 bgcolor=\'#336699\'><tr><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Name</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>PSSM</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Consensus</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>P-Value with Reference</font></th></tr>\n')
            clusterHtmlEnd = ''
            for member in bin1:
                # Plot the maxHit PSSM
                clusterHtml.write('<tr>')
                clusterHtml.write('<td bgcolor=\'#FFFFCC\'>'+member+'</td>') # Name
                pssmImgName = 'imgs/'+str(member)+'.png'
                #plotPssm(pssms[member], csDir+'/'+pssmImgName)
                clusterHtml.write('<td bgcolor=\'#FFFFCC\'><img src=\''+pssmImgName+'\'></td>') # PSSM Plot
                clusterHtml.write('<td bgcolor=\'#FFFFCC\'>'+pssms[member].getConsensusMotif()+'</td>') # Consensus
                clusterHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+str(binSimP[member])+'</center></td>')
                clusterHtml.write('</tr>\n')
                # Added so we can send the best hit PSSM to STAMP
                clusterHtmlEnd += '<div id=\''+str(member)+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
                clusterHtmlEnd += '<span class=\'gaggle-name hidden\'>'+str(member)+'</span>\n'
                clusterHtmlEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
                clusterHtmlEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[member].getMatrix()))+'x4</span>\n'
                clusterHtmlEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
                clusterHtmlEnd += '\n'.join(convertPSSM(pssms[member].getMatrix()))
                clusterHtmlEnd += '</div>\n</div>\n'
            clusterHtml.write('</table></center>\n'+clusterHtmlEnd+'</body>\n</html>')
            clusterHtml.close()
        for pssm1 in singletons:
            if float(pssms[pssm1].getEValue()) <= eValueThreshold:
                #outFile.write('\n'+str(pssm1)+','+pssms[pssm1].getConsensusMotif()+','+str(1)+','+'singletons')
                clusteringHtml.write('<tr><td bgcolor=\'#FFFFCC\'>'+pssm1+'</td>') # Singleton
                pssmImgName = 'imgs/'+str(pssm1)+'.png'
                plotPssm(pssms[pssm1], csDir+'/'+pssmImgName)
                clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><img src=\''+pssmImgName+'\'></td>') # Best hit PSSM Plot
                clusteringHtml.write('<td bgcolor=\'#FFFFCC\'>'+pssms[pssm1].getConsensusMotif()+'</td>') # Best hit Consensus
                clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center><b>Singleton</b></center></td></tr>') # Link to see all motifs in cluster
                # Added so we can send the best hit PSSM to STAMP
                clusteringHtmlEnd += '<div id=\''+str(pssm1)+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
                clusteringHtmlEnd += '<span class=\'gaggle-name hidden\'>'+str(pssm1)+'</span>\n'
                clusteringHtmlEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
                clusteringHtmlEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[pssm1].getMatrix()))+'x4</span>\n'
                clusteringHtmlEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
                clusteringHtmlEnd += '\n'.join(convertPSSM(pssms[pssm1].getMatrix()))
                clusteringHtmlEnd += '</div>\n</div>\n'
        #outFile.close()
        clusteringHtml.write('</table></center>\n'+clusteringHtmlEnd+'</body>\n</html>')
        print 'Done.\n'
    overallHtml.write('</table></center>\n')
    # Plot the barplot for the number of E-Values in each clustersize
    lp1 = plot(array(clusterSizes),array(goodEValsAll),'-ro')
    title('Cluster Size vs. E-Values <= '+str(eValueThreshold))
    xlabel('Cluster Size')
    ylabel('E-Values <= '+str(eValueThreshold))
    savefig('results/'+region+'/imgs/barPlot.png',format='png')
    close()

    # Make directories for saving output
    csDir = 'results/'+region+'/superBin'
    if not os.path.exists(csDir):
        os.mkdir(csDir)
    if not os.path.exists(csDir+'/imgs'):
        os.mkdir(csDir+'/imgs')
    
    # Make a pssms dictionary for TomTom runs
    pssms = {}
    motNs = {}
    for size in motifBins:
        for motif in motifBins[size]:
            pssms[str(size)+'-'+str(motif)] = motifBins[size][motif]['pssm']
            pssms[str(size)+'-'+str(motif)].setName(str(size)+'-'+str(motif))
            motNs[str(size)+'-'+str(motif)] = motifBins[size][motif]['numHits']
    
    # Do TOMTOM for bin vs bin from different cluster sizes
    # First make all files where each query file is one pssm and the target file is all remaining pssms, do the upper triangle
    pssmsKeys = pssms.keys()
    print 'Making Files ('+str(len(pssmsKeys))+')...'
    pssmNumDict = dict(zip(range(len(pssmsKeys)),pssmsKeys))
    for i in range(len(pssms.keys())):
        queryPssms = [deepcopy(pssms[pssmsKeys[i]])]
        queryPssms[0].setName(str(i))
        targetPssms = []
        pssmNames = deepcopy(pssmsKeys)
        pssmNames.remove(pssmsKeys[i])
        for other in pssmNames:
            targetPssms.append(pssms[other])
        makeFiles(nucFreqs=nucFreqs, queryPssms=queryPssms, targetPssms=targetPssms, num=i)
    print 'Done.'
    # Run this using all cores available
    cpus = cpu_count()
    print 'There are', cpus,'CPUs avialable.' 
    pool = Pool(processes=cpus)
    pool.map(runTomtom,range(len(pssms.keys())))
    print 'Done with Tomtom runs.\n'
    tomtomPValues = {}
    print 'Reading in Tomtom run...'
    for run in range(len(pssms.keys())):
        outputFile = open('tmp/tomtom_out/tomtom'+str(run)+'.out','r')
        output = outputFile.readlines()
        outputFile.close()
        # Now iterate through output and save data
        output.pop(0) # Get rid of header
        while len(output)>0:
            outputLine = output.pop(0).strip().split('\t')
            if len(outputLine)==9:
                id1 = pssmNumDict[int(outputLine[0])]
                id2 = outputLine[1]
                if not id1 in tomtomPValues:
                    tomtomPValues[id1] = { id2: { 'pValue': float(outputLine[3]), 'qValue':outputLine[4], 'offset':outputLine[2], 'overlap': outputLine[5], 'orientation': outputLine[8] } }
                else:
                    tomtomPValues[id1][id2] = { 'pValue': float(outputLine[3]), 'qValue':outputLine[4], 'offset':outputLine[2], 'overlap': outputLine[5], 'orientation': outputLine[8] }
                if id2 in tomtomPValues and id1 in tomtomPValues[id2]:
                    if tomtomPValues[id1][id2] > tomtomPValues[id2][id1]:
                        tomtomPValues[id1][id2] = tomtomPValues[id2][id1]
                    else:
                        tomtomPValues[id2][id1] = tomtomPValues[id1][id2]
    print 'Done.\n'
    
    ## Agglomerative complete linkage clustering algorithm
    #   Step 1. Find minimum pairwise PSSM p-value between all clusters
    #   Step 2. Combine clusters
    #   Step 3. Continue until p-value between clusters is greater than pValueThreshold
    print 'SuperBin: Agglomerative complete linkage clustering start...'
    clusters = [[i] for i in deepcopy(pssms.keys())]
    done = 0
    while 1:
        # Find minimum pairwise PSSM p-value (complete linkage)
        clusters2 = deepcopy(clusters)
        minP = ['', '', float(1)]
        for c1 in clusters:
            clusters2.pop(0)
            for c2 in clusters2:
                testMe = completeLinkage(c1, c2, tomtomPValues, pValueThreshold)
                if testMe<minP[2]:
                    minP = [c1, c2, testMe]
        if minP[2]<pValueThreshold:
            clusters.remove(minP[0])
            clusters.remove(minP[1])
            clusters.append(minP[0]+minP[1])
        else:
            break
    print len(clusters), clusters
    print 'Done.\n'

    # Remove singletons and sort bins
    print 'SuperBin: Sort bins and remove singletons...'
    binned = []
    binSimP = {}
    singletons = []
    for bin1 in clusters:
        if len(bin1)>1:
            # Pick the exemplar motif
            if len(bin1)>2:
                # Pick the most representative PSSM
                pssmHits = {}
                for pssm1 in bin1:
                    pssmHits[pssm1] = 0
                    for pssm2 in bin1:
                        if not pssm1==pssm2 and tomtomPValues[pssm1][pssm2]['pValue'] <= pValueThreshold:
                            pssmHits[pssm1] += 1
                maxHit = getMaxHit(pssmHits)
                tmpBin = deepcopy(bin1)
                tmpBin.remove(maxHit)
                # Sort the remaining based on similarity to maximum hit
                basedOn = []
                for i in tmpBin:
                    basedOn.append(tomtomPValues[maxHit][i]['pValue'])
                sortedOut = qsortBasedOn(tmpBin,basedOn)
                binned.append(deepcopy([maxHit]+sortedOut[0]))
                binSimP.update(dict(zip([maxHit]+sortedOut[0],['Reference']+sortedOut[1])))
            else:
                binned.append(deepcopy(bin1))
                binSimP.update(dict(zip(bin1,['-','-'])))
        else:
            singletons.append(bin1[0])
    binned = qsortBins(binned)
    print 'Done.\n'

    # Write out the superBins
    print 'Write out SuperBins...'
    clusteringHtml = open(csDir+'/clustering.html','w')
    clusteringHtml.write('<html>\n<head><title>SuperBin MEME Results</title>\n<style type=\'text/css\'>\n.hidden {\ndisplay: none;\n }\n</style>\n</head>\n<body style=\'font-family: arial, verdana, sans-serif\'>\n')
    clusteringHtml.write('<center><table border=1 bgcolor=\'#336699\'><tr><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Best SuperMatch</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>PSSM</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Consensus</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Number of SuperMatches</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Overall Motif Matches</font></th></tr>\n')
    clusteringHtmlEnd = ''
    for bin1 in binned:
        clusteringHtml.write('<tr>')
        # Plot the maxHit PSSM
        clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+bin1[0]+'</center></td>') # Best hit
        pssmImgName = 'imgs/'+str(bin1[0])+'.png'
        plotPssm(pssms[bin1[0]], csDir+'/'+pssmImgName)
        clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><img src=\''+pssmImgName+'\'></td>') # Best hit PSSM Plot
        clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+pssms[bin1[0]].getConsensusMotif()+'</center></td>') # Best hit Consensus
        clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center><b><a href=\'bin'+bin1[0]+'.html\'>'+str(len(bin1))+'</a></b></center></td>') # Link to see all motifs in cluster
        totalN = 0
        for mot in bin1:
            totalN += motNs[mot]
        clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><center>'+str(totalN)+'</center></td>') # Link to see all motifs in cluster
        clusteringHtml.write('</tr>\n')
        # Added so we can send the best hit PSSM to STAMP
        clusteringHtmlEnd += '<div id=\''+str(bin1[0])+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
        clusteringHtmlEnd += '<span class=\'gaggle-name hidden\'>'+str(bin1[0])+'</span>\n'
        clusteringHtmlEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
        clusteringHtmlEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[bin1[0]].getMatrix()))+'x4</span>\n'
        clusteringHtmlEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
        clusteringHtmlEnd += '\n'.join(convertPSSM(pssms[bin1[0]].getMatrix()))
        clusteringHtmlEnd += '</div>\n</div>\n'
        # Now build each clusters html file
        clusterHtml = open(csDir+'/bin'+bin1[0]+'.html','w')
        clusterHtml.write('<html>\n<head><title>SuperBin '+bin1[0]+' of MEME Results for Cluster Size '+str(clusterSize)+'</title>\n<style type=\'text/css\'>\n.hidden {\ndisplay: none;\n }\n</style>\n</head>\n<body style=\'font-family: arial, verdana, sans-serif\'>\n')
        clusterHtml.write('<center><table border=1 bgcolor=\'#336699\'><tr><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Name</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>PSSM</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>Consensus</font></th><th bgcolor=\'#336699\'><font color=\'#FFFFCC\'>P-Value with Reference</font></th></tr>\n')
        clusterHtmlEnd = ''
        for member in bin1:
            # Plot the maxHit PSSM
            clusterHtml.write('<tr>')
            clusterHtml.write('<td bgcolor=\'#FFFFCC\'>'+member+'</td>') # Name
            pssmImgName = 'imgs/'+str(member)+'.png'
            plotPssm(pssms[member], csDir+'/'+pssmImgName)
            clusterHtml.write('<td bgcolor=\'#FFFFCC\'><img src=\''+pssmImgName+'\'></td>') # PSSM Plot
            clusterHtml.write('<td bgcolor=\'#FFFFCC\'>'+pssms[member].getConsensusMotif()+'</td>') # Consensus
            clusterHtml.write('<td bgcolor=\'#FFFFCC\'>'+str(binSimP[member])+'</td>') # Consensus
            clusterHtml.write('</tr>\n')
            # Added so we can send the best hit PSSM to STAMP
            clusterHtmlEnd += '<div id=\''+str(member)+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
            clusterHtmlEnd += '<span class=\'gaggle-name hidden\'>'+str(member)+'</span>\n'
            clusterHtmlEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
            clusterHtmlEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[member].getMatrix()))+'x4</span>\n'
            clusterHtmlEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
            clusterHtmlEnd += '\n'.join(convertPSSM(pssms[member].getMatrix()))
            clusterHtmlEnd += '</div>\n</div>\n'
        clusterHtml.write('</table></center>\n'+clusterHtmlEnd+'</body>\n</html>')
        clusterHtml.close()
    for pssm1 in singletons:
        if float(pssms[pssm1].getEValue()) <= eValueThreshold:
            clusteringHtml.write('<tr><td bgcolor=\'#FFFFCC\'>'+pssm1+'</td>') # Singleton
            pssmImgName = 'imgs/'+str(pssm1)+'.png'
            plotPssm(pssms[pssm1], csDir+'/'+pssmImgName)
            clusteringHtml.write('<td bgcolor=\'#FFFFCC\'><img src=\''+pssmImgName+'\'></td>') # Best hit PSSM Plot
            clusteringHtml.write('<td bgcolor=\'#FFFFCC\'>'+pssms[pssm1].getConsensusMotif()+'</td>') # Best hit Consensus
            clusteringHtml.write('<td colspan=2 bgcolor=\'#FFFFCC\'><center><b>Singleton</b></center></td></tr>') # Link to see all motifs in cluster
            # Added so we can send the best hit PSSM to STAMP
            clusteringHtmlEnd += '<div id=\''+str(pssm1)+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
            clusteringHtmlEnd += '<span class=\'gaggle-name hidden\'>'+str(pssm1)+'</span>\n'
            clusteringHtmlEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
            clusteringHtmlEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[pssm1].getMatrix()))+'x4</span>\n'
            clusteringHtmlEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
            clusteringHtmlEnd += '\n'.join(convertPSSM(pssms[pssm1].getMatrix()))
            clusteringHtmlEnd += '</div>\n</div>\n'
    clusteringHtml.write('</table></center>\n'+clusteringHtmlEnd+'</body>\n</html>')
    print 'Done.\n'
overallHtml.write('</body>\n</html>')
overallHtml.close()
