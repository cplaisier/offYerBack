#################################################################
# @Program: tomtom.py                                           #
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

import os
from subprocess import *

# A class designed to compute the difference between two PSSM matrices.
#
# Variables:
# 
# 
# Functions:
# getScore(self,queryPssm,targetPssm)
# getAllScores(self)
#
class tomtom:
    # Initialize and start the run
    def __init__(self,queryPssms,targetPssms,nucFreqs,strands='+ -',distMeth='ed',qThresh=1,minOverlap=0):
        # Set some parameters for analysis
        self.queryPssms = queryPssms
        self.targetPssms = targetPssms
        self.nucFreqs = nucFreqs
        self.strands = strands
        self.distMeth = distMeth
        self.qThresh = qThresh
        self.minOverlap = minOverlap
        # Setup for a tomtom run
        self.makeFiles()
        # Run tomtom and save data into local data structures
        print 'Starting TOMTOM...'
        self.runTomTom()
        print 'Done.\n'

    # Make the files for a TomTom run
    def makeFiles(self):
        if not os.path.exists('tmp'):
            os.mkdir('tmp')
        # Header crap
        memeHeader = ''
        memeHeader += 'MEME version 3.0\n\n'
        memeHeader += 'ALPHABET= ACGT\n\n'
        # Here is where we tell it what strand: for miRNAs this would just be '+'
        memeHeader += 'strands: '+self.strands+'\n\n'
        memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
        memeHeader += 'A '+str(round(float(self.nucFreqs['A']),3))+' C '+str(round(float(self.nucFreqs['C']),3))+' G '+str(round(float(self.nucFreqs['G']),3))+' T '+str(round(float(self.nucFreqs['T']),3))
        #extras = ['R','Y','K','M','S','W','B','D','H','V','N']
        #for e in extras:
        #    memeHeader += ' '+str(e)+' 0.000'
        # Make query PSSM file
        queryFile = open('tmp/query.tomtom','w')
        queryFile.write(memeHeader)
        for pssm in self.queryPssms:
            queryFile.write('\n\n'+pssm.getMemeFormatted())
        queryFile.close()
        # Make target PSSM file
        targetFile = open('tmp/target.tomtom','w')
        targetFile.write(memeHeader)
        for pssm in self.targetPssms:
            targetFile.write('\n\n'+pssm.getMemeFormatted())
        targetFile.close()

    # Run TomTom on the files
    def runTomTom(self):
        if not os.path.exists('tmp/tomtom_out'):
            os.mkdir('tmp/tomtom_out')
        # Arguments for tomtom
        tomtomArgs = ' -query tmp/query.tomtom -target tmp/target.tomtom -dist '+str(self.distMeth)+' -o tmp/tomtom_out -text -q-thresh '+str(self.qThresh)+' -min-overlap '+str(self.minOverlap)+' -verbosity 0'
        print tomtomArgs
        #p = Popen("tomtom" + tomtomArgs, shell=True)
        #sts = os.waitpid(p.pid, 0)
        errOut = open('tmp/stderr.out','w')
        tomtomProc = Popen("tomtom" + tomtomArgs, shell=True,stdout=PIPE,stderr=errOut)
        output = tomtomProc.communicate()[0].split('\n')
        errOut.close()
        # Now iterate through output and save data
        output.pop(0) # Get rid of header
        self.scoreMatrix = {}
        while len(output)>0:
            outputLine = output.pop(0).strip().split('\t')
            if len(outputLine)==9:
                if not outputLine[0] in self.scoreMatrix:
                    self.scoreMatrix[outputLine[0]] = { outputLine[1]: { 'optimalOffset':outputLine[2], 'pValue':outputLine[3], 'qValue':outputLine[4], 'overlap':outputLine[5], 'queryConsensus':outputLine[6], 'targetConsensus':outputLine[7], 'orientation':outputLine[8] } }
                else:
                    self.scoreMatrix[outputLine[0]][outputLine[1]] = { 'optimalOffset':outputLine[2], 'pValue':outputLine[3], 'qValue':outputLine[4], 'overlap':outputLine[5], 'queryConsensus':outputLine[6], 'targetConsensus':outputLine[7], 'orientation':outputLine[8] }
    
    # Get all pair-wise scores
    def getAllScores(self):
        return self.scoreMatrix

    # Get a specific pair-wise ananlysis result 
    def getScore(self,queryPssmName,targetPssmName):
        scores = self.getAllScores()
        if not queryPssmName in scores:
            return ['No such query pssm name!']
        else:
            if not targetPssmName in scores[queryPssmName]:
                return ['No such target!']
            else:
                return scores[queryPssmName][targetPssmName]

