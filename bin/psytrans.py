#!/usr/bin/env python

# Script was modified from original to use python3

import argparse
import array
import gzip
import logging
import math
import os
import os.path
import random
import shutil
import subprocess
import string
import sys
import tempfile
import threading
import traceback

if(sys.hexversion < 0x03000000):
    import Queue
else:
    import queue as Queue


########################
########################
### Global constants ###
########################
########################

PSYTRANS_VERSION = '1.0.0'

HOST_NAME     = 'species1'
SYMB_NAME     = 'species2'
DB_NAME       = 'psytrans_refenrece_db'
DB_FASTA      = DB_NAME + '.fasta'
BLAST_FILE    = HOST_NAME + SYMB_NAME + '_blastResults.txt'
BLAST_SORT    = HOST_NAME + SYMB_NAME + '_blastClassification.txt'
HOST_CODE     = 1
SYMB_CODE     = 2
HOST_CODE_STR = str(HOST_CODE)
SYMB_CODE_STR = str(SYMB_CODE)
UNKNOWN_CODE  = 3

#Blast classification variables
CLASSIFICATION_MIN_BIT_RATIO = 2
CLASSIFICATION_MIN_BIT_DELTA = 100
TRAINING_MIN_BIT_RATIO       = 5
TRAINING_MIN_BIT_DELTA       = 400

# Proportion of sequences to keep for training if less than args.numberOfSeq
# sequences can be used for training
TRAINING_PROPORTION = 0.75

# SVM global variables
SVM_FOLD   = 5
SVM_CSTART = -5
SVM_CEND   = 15
SVM_CSTEP  = 2
SVM_GSTART = 3
SVM_GEND   = -15
SVM_GSTEP  = -2

HOST_TRAINING = HOST_NAME + '_training.fasta'
HOST_TESTING  = HOST_NAME + '_testing.fasta'
SYMB_TRAINING = SYMB_NAME + '_training.fasta'
SYMB_TESTING  = SYMB_NAME + '_testing.fasta'
BINARIES_DIR  = 'binaries'

LETTERS = ('A', 'T', 'G', 'C')

COMPLEMENT_TABLE = str.maketrans('ATGCatgc', 'TACGtacg')

####################################################################
####################################################################
### Class to store the various paths used throughout the program ###
####################################################################
####################################################################

class PsyTransOptions:
    """This class consists of attributes to allow database and file paths to be
    obtained conveniently and consistently."""

    def __init__(self, args):
        self.args              = args
        self.dbPath            = None
        self.fastaDbPath       = None
        self.blastResultsPath  = None
        self.suffix            = None
        self.logName           = None
        self.inputFile         = None
        self.trainFile         = None
        self.testFile          = None
        self.hostTrainPath     = None
        self.hostTestPath      = None
        self.symbTrainPath     = None
        self.symbTestPath      = None
        self.blastSortPath     = None
        self.SVMOutPath        = None
        self.tempDir           = None
        self.chunkList         = []
        self.threadBlastList   = []
        self.createTempDir()

    def createTempDir(self):
        self.tempDir = tempfile.mkdtemp(prefix='psytrans_',
                                        suffix='_temp',
                                        dir=self.args.tempDir)

    def getDbPath(self):
        """Return the path to the blast database"""
        if not self.dbPath:
            dbPath      = os.path.join(self.tempDir, DB_NAME)
            self.dbPath = os.path.abspath(dbPath)
        return self.dbPath

    def getFastaDbPath(self):
        """Return the path of the fasta file for the blast database"""
        if not self.fastaDbPath:
            self.fastaDbPath = os.path.join(self.tempDir, DB_FASTA)
        return self.fastaDbPath

    def getChunkList(self):
        """Return the list of path to the fasta chunks"""
        if not self.chunkList:
            blastInput  = self.args.queries
            fastaPrefix = os.path.basename(blastInput)
            for i in range(self.args.nbThreads):
                fastaChunkName = '%s_chunk_%06d' % (fastaPrefix, i)
                fastaChunkName = os.path.join(self.tempDir, fastaChunkName)
                self.chunkList.append(fastaChunkName)
        return self.chunkList

    def getThreadBlastList(self):
        """Return the list of output paths of the multi threaded blast"""
        if not self.threadBlastList:
            for i in range(self.args.nbThreads):
                blastThreadFile = '%s.%06d' % (BLAST_FILE, i)
                blastThreadPath = os.path.join(self.tempDir, blastThreadFile)
                self.threadBlastList.append(blastThreadPath)
        return self.threadBlastList

    def getBlastResultsPath(self):
        """Return the path to the blast results"""
        if self.args.blastResults:
            return self.args.blastResults
        if not self.blastResultsPath:
            self.blastResultsPath = os.path.join(self.tempDir, BLAST_FILE)
        return self.blastResultsPath

    def getCheckPointPath(self, dFile):
        """Return the path for the chekpoint (.done) file"""
        return os.path.join(self.tempDir, dFile)

    def createCheckPoint(self, cpFile):
        """Create the checkpoint file"""
        path = self.getCheckPointPath(cpFile)
        open(path, "w").close()

    def checkPoint(self, dFile):
        """Check if a particular checkpoint has been created"""
        path = self.getCheckPointPath(dFile)
        return os.path.exists(path)

    def _getNumberOfSequences(self):
        """Return the length part of the kmer file name"""
        if self.args.numberOfSeq == 0:
            length = 'all'
        else:
            length = self.args.numberOfSeq
        return str(length)

    def _getSuffix(self):
        """Create the suffix of the SVM input files"""
        if not self.suffix:
            suffix = self._getNumberOfSequences()
            self.mink = str(self.args.minWordSize)
            self.maxk = str(self.args.maxWordSize)
            self.suffix = suffix + '_c' + self.mink + '_k' + self.maxk
        return self.suffix

    def getLogPath(self):
        if not self.logName:
            self.logName = 'psytrans_' + self._getSuffix() + '.log'
        return self.logName

    def getTrainPath(self):
        """Return the path of the SVM training file"""
        if not self.trainFile:
            fName = self._getSuffix()
            self.trainFile = 'Training' + '_' + fName + '.txt'
        return self.trainFile

    def getTestPath(self):
        """Return the path of the SVM testing file"""
        if not self.testFile:
            fName = self._getSuffix()
            self.testFile = 'Testing' + '_' + fName + '.txt'
        return str(self.testFile)

    def getHostTrainPath(self):
        """Return the path of the host training sequences"""
        if not self.hostTrainPath:
            self.hostTrainPath = os.path.join(self.tempDir, HOST_TRAINING)
        return self.hostTrainPath

    def getHostTestPath(self):
        """Return the path of the host testing sequences"""
        if not self.hostTestPath:
            self.hostTestPath = os.path.join(self.tempDir, HOST_TESTING)
        return self.hostTestPath

    def getSymbTrainPath(self):
        """Return the path of the symbiont training sequences"""
        if not self.symbTrainPath:
            self.symbTrainPath = os.path.join(self.tempDir, SYMB_TRAINING)
        return self.symbTrainPath

    def getSymbTestPath(self):
        """Return the path of the symbiont testing sequences"""
        if not self.symbTestPath:
            self.symbTestPath = os.path.join(self.tempDir, SYMB_TESTING)
        return self.symbTestPath

    def getBlastSortPath(self):
        """Get the path of the sequences sorted using blast"""
        if not self.blastSortPath:
            self.blastSortPath = os.path.join(self.tempDir, BLAST_SORT)
        return self.blastSortPath

    def getSVMOutPath(self):
        """Get the SVM output path"""
        if not self.SVMOutPath:
            fName = self._getSuffix()
            self.SVMOutPath = fName + '.out'
            self.SVMOutPath = os.path.join(self.tempDir, self.SVMOutPath)
        return self.SVMOutPath

######################
######################
### Misc utilities ###
######################
######################

def iterFasta(path):
    """Iterates over the sequences of a fasta file"""

    logging.info("Loading fasta files from %s" % path)
    name = None
    seq = []
    if path.endswith('.gz') or path.endswith('.gz"'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name = line[1:]
            seq  = []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))
    handle.close()

def seqCount(path):
    """Counts the number of sequences in a fasta file"""

    c = 0
    if path.endswith('.gz'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        if line.startswith(">"):
            c += 1
    return c

def revComp(seq):
    return seq.translate(COMPLEMENT_TABLE)[::-1]

#####################################
#####################################
### Make training set using BLAST ###
#####################################
#####################################

def writeDatabase(options, fastaPath):
    """Write Host and Symbiont sequences with standardised names to a new fasta
    file"""

    logging.info('Creating Database.')
    targetPath = open(fastaPath, "w")
    #Writing Host Sequences to target database
    i = 0
    for name, seq in iterFasta(options.args.species1):
        i += 1
        targetPath.write('>%s_%d\n%s\n' % (HOST_NAME, i, seq))
    i = 0
    for name, seq in iterFasta(options.args.species2):
        i += 1
        targetPath.write('>%s_%d\n%s\n' % (SYMB_NAME, i, seq))
    targetPath.close()
    options.createCheckPoint('writedatabase.done')

def makeDB(options):
    """Build the blast database in the temporary folder"""

    dbPath      = options.getDbPath()
    fastaPath   = options.getFastaDbPath()
    logPath     = dbPath + '.log'
    makeblastdb = checkExecutable('makeblastdb')
    dbType      = 'nucl'
    if options.args.blastType == 'blastx':
        dbType  = 'prot'
    makeDBCmd   = [makeblastdb,
                   '-title',
                   DB_NAME,
                   '-in',
                   fastaPath,
                   '-dbtype',
                   dbType,
                   '-out ',
                   dbPath,
                   '-logfile',
                   logPath]
    makeDBCmd   = ' '.join(makeDBCmd)
    submakeDB   = subprocess.call(makeDBCmd, shell=True)
    if not submakeDB == 0:
        logging.error('[ERROR] Failed to create blast database')
        sys.exit(1)
    options.createCheckPoint('makeDB.done')

def splitBlastInput(options):
    """Split the input fasta file into chunks for parallel blast searches"""

    logging.info('Splitting sequences into %d chunks' % options.args.nbThreads)
    chunkList = options.getChunkList()
    handles   = []
    for i in range(options.args.nbThreads):
        handle = open(chunkList[i], "w")
        handles.append(handle)
    #writing to each chunk .fasta
    i = 0
    for name, seq in iterFasta(options.args.queries):
        handles[i % options.args.nbThreads].write('>%s\n%s\n' % (name, seq))
        i += 1
    for i in range(options.args.nbThreads):
        handles[i].close()

def runBlast(options, threadId):
    """Invoke the blast command. The output format of the result by default is
    set to '6' (tab-seaparated without comments)."""

    #Define BLAST variables
    logging.info('Performing Blast search with thread %d' % threadId)
    eVal         = '%.2e' % options.args.maxBestEvalue
    dbPath       = options.getDbPath()
    blastOutput  = options.getThreadBlastList()[threadId]
    blast        = checkExecutable(options.args.blastType)
    blastCmd     = [blast,
                    '-evalue',
                    eVal,
                    '-query',
                    options.getChunkList()[threadId],
                    '-db',
                    dbPath,
                    '-outfmt 6',
                    '-out',
                    blastOutput]
    blastCmd     = ' '.join(blastCmd)
    retCode      = subprocess.call(blastCmd, shell=True)
    if not retCode == 0:
        logging.error('[ERROR] Failed to excecute blast command')
        sys.exit(1)

def mergeBlastOutput(options):
    """Merge the output from the blast searches"""

    logging.info('Merging Blast results')
    blastOut       = options.getBlastResultsPath()
    blastOutHandle = open(blastOut, "w")
    for i in range(options.args.nbThreads):
        threadPath   = options.getThreadBlastList()[i]
        threadHandle = open(threadPath)
        for line in threadHandle:
            blastOutHandle.write(line)
        threadHandle.close()
    blastOutHandle.close()

def runBlastThreads(options):
    """Split the queries into chunks, run blast and merge the results"""

    logging.info('Launching threaded Blast search')
    splitBlastInput(options)
    threads = []
    for i in range(options.args.nbThreads):
        t = threading.Thread(target=runBlast, args=(options, i))
        threads.append(t)
        t.start()
    for i in range(options.args.nbThreads):
        threads[i].join()
    mergeBlastOutput(options)
    options.createCheckPoint('runBlast.done')

def processBlastResults(options):
    """Wraps the whole blast processing in a single function to make checkpointing easier"""

    querries = parseBlast(options)
    trainingClassification, blastClassification = classifyFromBlast(options, querries)
    seqSplit(options, trainingClassification)
    writeBlastClassifications(options, blastClassification)
    options.createCheckPoint('parseBlast.done')

def parseBlast(options):
    """Parse the blast results to be used later to prepare the training and
    testing set with unambiguously classified sequences.
    The expected format is tabular (-outfmt 6 or -outfmt 7)."""

    logging.info('Parsing blast results')
    path = options.getBlastResultsPath()
    if path.endswith('.gz'):
        handle = gzip.open(path)
    else:
        handle = open(path, "r")
    querries = {}
    n        = 0
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '#':
            continue
        fields = line.split()
        qName  = fields[0]
        hName  = fields[1]
        evalue = float(fields[10])
        bitscore = float(fields[11])
        if not qName in querries:
            querries[qName] = []
            n += 1

        hit = (hName, evalue, bitscore)
        querries[qName].append(hit)
    logging.info('Parsed %d blast records' % n)
    logging.info('Found %d queries hits' % len(querries))
    handle.close()
    return querries

def classifyFromBlast(options, querries):
    """Classify blast results into ambiguous and unambiguous sequences from
    host and symbiont"""

    def sortHits(h1, h2):
        """Sorts hits by evalue ( the second field in this 3-tuple)"""

        if h1[1] > h2[1]:
            return 1
        if h1[1] < h2[1]:
            return -1
        return 0

    logging.info('Classifying using Blast results')
    trainingClassification = {}
    blastClassification    = {}
    hostTrained    = 0
    symbTrained    = 0
    hostClassified = 0
    symbClassified = 0
    for qName in querries:
        hits = querries[qName]
        hits.sort(key= lambda x: sortHits(x,x))
        hasCoral        = False
        hasZoox         = False
        coralBestEvalue = -1
        zooxBestEvalue  = -1
        coralBestScore  = 0
        zooxBestScore   = 0
        for hName, evalue, bitscore in hits:
            if hName.startswith(HOST_NAME):
                if not hasCoral:
                    hasCoral        = True
                    coralBestEvalue = evalue
                    coralBestScore  = bitscore
            else :
                if not hasZoox:
                    hasZoox        = True
                    zooxBestEvalue = evalue
                    zooxBestScore  = bitscore
        if hasCoral and not hasZoox and coralBestEvalue <= options.args.maxBestEvalue:
            trainingClassification[qName] = HOST_CODE
            blastClassification[qName]    = HOST_CODE
            hostTrained                  += 1
            hostClassified               += 1
        elif hasZoox and not hasCoral and zooxBestEvalue <= options.args.maxBestEvalue:
            trainingClassification[qName] = SYMB_CODE
            blastClassification[qName]    = SYMB_CODE
            symbTrained                  += 1
            symbClassified               += 1
        if hasZoox and hasCoral:
            # TODO these two conditionals could be factorised into one
            if coralBestScore > zooxBestScore:
                scoreRatio = float(coralBestScore) / float(zooxBestScore)
                scoreDelta = coralBestScore - zooxBestScore
                if scoreRatio > CLASSIFICATION_MIN_BIT_RATIO and scoreDelta > CLASSIFICATION_MIN_BIT_DELTA:
                    blastClassification[qName] = HOST_CODE
                    hostClassified            += 1
                if scoreRatio > TRAINING_MIN_BIT_RATIO and scoreDelta > TRAINING_MIN_BIT_DELTA:
                    trainingClassification[qName] = HOST_CODE
                    hostTrained                  += 1
            elif coralBestScore < zooxBestScore:
                scoreRatio = float(zooxBestScore) / float(coralBestScore)
                scoreDelta = zooxBestScore - coralBestScore
                if scoreRatio > CLASSIFICATION_MIN_BIT_RATIO and scoreDelta > CLASSIFICATION_MIN_BIT_DELTA:
                    blastClassification[qName] = SYMB_CODE
                    symbClassified            += 1
                if scoreRatio > TRAINING_MIN_BIT_RATIO and scoreDelta > TRAINING_MIN_BIT_DELTA:
                    trainingClassification[qName] = SYMB_CODE
                    symbTrained                  += 1
    logging.info('Found %d unambiguous hits' % len(trainingClassification))
    logging.info('Found %d host only hits' % hostTrained)
    logging.info('Found %d symbiont only hits' % symbTrained)
    logging.info('Found %d likely host hits' % hostClassified)
    logging.info('Found %d likely symbiont hits' % symbClassified)
    return trainingClassification, blastClassification

def seqSplit(options, trainingClassification):
    """Write the unambiguously classified sequences into four fasta files:
    training.fasta for host sequences, testing.fasta for host sequences,
    training.fasta for symb sequences and testing.fasta for symb sequences."""

    logging.info('Splitting training and testing sequences')
    tooSmall                 = 0
    longSeqs                 = {}
    longSeqs[SYMB_CODE]      = []
    longSeqs[HOST_CODE]      = []
    outHandles               = {}
    outHandles[HOST_CODE]    = (open(options.getHostTrainPath(), "w"),
                                open(options.getHostTestPath(), "w"))
    outHandles[SYMB_CODE]    = (open(options.getSymbTrainPath(), "w"),
                                open(options.getSymbTestPath(), "w"))

    for name, seq in iterFasta(options.args.queries):
        size = len(seq)
        if size < options.args.minSeqSize:
            tooSmall += 1
            continue
        name = name.split()[0]
        if name in trainingClassification:
            classification = trainingClassification[name]
            longSeqs[classification].append((name, seq))

    # Shuffle so that the results are independent of the order of the input sequences
    random.shuffle(longSeqs[SYMB_CODE])
    random.shuffle(longSeqs[HOST_CODE])

    for classification in (HOST_CODE, SYMB_CODE):
        seqs           = longSeqs[classification]
        nSeqs          = len(seqs)
        trainingMaxIdx = options.args.numberOfSeq
        if nSeqs < 2 * options.args.numberOfSeq:
            trainingMaxIdx = int(nSeqs * TRAINING_PROPORTION)
        # Training file
        for i in range(trainingMaxIdx):
            outHandles[classification][0].write('>%s\n%s\n' % seqs[i])
        # Testing file
        for i in range(trainingMaxIdx, 2 * options.args.numberOfSeq):
            outHandles[classification][1].write('>%s\n%s\n' % seqs[i])

    for classification in (HOST_CODE, SYMB_CODE):
        outHandles[classification][0].close()
        outHandles[classification][1].close()
    logging.info('When creating training and testing sets, %d sequences were smaller than %d and ignored' % \
                 (tooSmall, options.args.minSeqSize))
    # FIXME log the number of sequences writen for training and testing of each species

def writeBlastClassifications(options, blastClassification):
    """Write a file with the blast classification"""

    handle = open(options.getBlastSortPath(), "w")
    for name in blastClassification:
        classification = blastClassification[name]
        handle.write('%s\t%d\n' % (name, classification))

############################
############################
### Compute Kmer vectors ###
############################
############################

def prepareMaps(k, maxk, kmers):
    """Prepares the kmer maps for the specified kmer range"""

    if k == maxk:
        n        = 0
        kmer2int = {}
        for kmer in kmers:
            kmer2int[kmer] = n
            n += 1
        return kmer2int
    newKmers = []
    for kmer in kmers:
        for letter in LETTERS:
            newKmers.append(kmer + letter)
    kmers = newKmers
    return prepareMaps(k + 1, maxk, kmers)

def computeKmers(options, path, outfile, code, mode):
    """Compute the kmer counts throughout the kmer range for each sequence, and
    write the output to a file.  Each kmer counts will be scaled accordingly
    with the sequence size."""

    logging.info('Computing kmers for %s' % path)
    # Prepare all maps
    kMin    = options.args.minWordSize
    kMax    = options.args.maxWordSize
    maps    = []
    logging.info('Preparing kmer maps')
    for i in range(kMin, kMax + 1):
        maps.append(prepareMaps(0, i, ['']))
    # Initialise output
    out     = outfile
    outPath = os.path.join(options.tempDir, out)
    handle  = open(outPath, mode)
    # Initialise counts
    counts  = {}
    for i in range(kMin, kMax + 1):
        counts[i] = array.array('d', [0 for x in range(4 ** i)])
    # Iterate over sequences
    nSeqs   = 0
    for name, seq in iterFasta(path):
        size   = len(seq)
        n      = 0
        handle.write('%d' % code)
        seqs = (seq, )
        if options.args.bothStrands:
            seqs = (seq, revComp(seq))
        # For each kmer value
        for i in range(kMin, kMax + 1):
            kCounts  = counts[i]
            # For each strand
            for seq in seqs:
                # For each word in the sequence
                for j in range(size - i + 1):
                    word = seq[j:j + i]
                    kMap = maps[i - kMin]
                    idx  = kMap.get(word, None)
                    # If the word contains characters other than ATGC
                    if idx is None:
                        continue
                    kCounts[idx] += 1
            kCountsSum = sum(kCounts)
            for j in range(len(kCounts)):
                kCounts[j] /= kCountsSum
            for j in kCounts:
                n += 1
                if j != 0:
                    handle.write(' %d:%.3e' % (n, j))
        handle.write('\n')
        nSeqs += 1
        # Reset counts
        for i in range(kMin, kMax + 1):
            for j in range(len(counts[i])):
                counts[i][j] = 0
    # Trace
    logging.info('Processed %d sequences' % nSeqs)
    handle.close()

def prepareTrainingKmers(options, kmerTrainPath, kmerTestPath):
    """Compute the kmer counts for the training and testing sequences.  The
    function outputs two files: a training file and a testing file to be used
    as inputs for the SVM training."""

    logging.info('Preparing kmers for training')
    hostTrainPath = options.getHostTrainPath()
    hostTestPath  = options.getHostTestPath()
    symbTrainPath = options.getSymbTrainPath()
    symbTestPath  = options.getSymbTestPath()
    computeKmers(options, hostTrainPath, kmerTrainPath, HOST_CODE, "w")
    computeKmers(options, hostTestPath, kmerTestPath, HOST_CODE, "w")
    computeKmers(options, symbTrainPath, kmerTrainPath, SYMB_CODE, "a")
    computeKmers(options, symbTestPath, kmerTestPath, SYMB_CODE, "a")
    options.createCheckPoint('kmers.done')

##################################################################
##################################################################
### SVM computations, based on svm-easy / svm-grid from libsvm ###
##################################################################
##################################################################

def doSVMEasy(options, kmerTrain, kmerTest):
    """Scale the input, optimise the SVM parameters and test the predictions.
    the This is roughly the equivalent of the svm-easy script from libsvm."""

    logging.info('Starting SVM training')
    kmerTrain       = os.path.join(options.tempDir, kmerTrain)
    kmerTest        = os.path.join(options.tempDir, kmerTest)
    svmTrain        = checkExecutable('svm-train')
    svmPredict      = checkExecutable('svm-predict')
    svmScale        = checkExecutable('svm-scale')
    scaledFile      = kmerTrain + '.scale'
    modelFile       = kmerTrain + '.model'
    rangeFile       = kmerTrain + '.range'
    scaledTestFile  = kmerTest  + '.scale'
    predictTestFile = kmerTest  + '.predict'
    resultLog       = kmerTrain + '_svmPredict.log'
    cmdScale        = [svmScale,
                       '-s',
                       rangeFile,
                       kmerTrain,
                       '>',
                       scaledFile]
    cmdScale        = ' '.join(cmdScale)
    subprocess.call(cmdScale, shell=True)
    c, g, rate      = doSVMGrid(options, scaledFile)
    cmdTrain        = [svmTrain,
                       '-c',
                       str(c),
                       '-g',
                       str(g),
                       scaledFile,
                       modelFile]
    cmdTrain        = ' '.join(cmdTrain)
    subprocess.call(cmdTrain, shell=True)
    cmdScale        = [svmScale,
                       '-r',
                       rangeFile,
                       kmerTest,
                       '>',
                       scaledTestFile]
    cmdScale        = ' '.join(cmdScale)
    subprocess.call(cmdScale, shell=True)
    cmdPredict      = [svmPredict,
                       scaledTestFile,
                       modelFile,
                       predictTestFile]
    cmdPredict      = ' '.join(cmdPredict)
    resultHandle    = open(resultLog, "w")
    subprocess.call(cmdPredict, shell=True, stdout=resultHandle)
    #Adding classification-result to logger
    for line in open(resultLog, "r"):
        if line.startswith('Accuracy'):
            logging.info(line.strip())
            break
    resultHandle.close()
    logging.info('Prediction in: %s' % predictTestFile)
    options.createCheckPoint('svm.done')

def calculateSVMGridJobs():
    """Calculate the coordinates of the search space"""

    def rangeF(begin, end, step):
        seq = []
        while True:
            if step > 0 and begin > end:
                break
            if step < 0 and begin < end:
                break
            seq.append(begin)
            begin = begin + step
        return seq

    def permuteSequence(seq):
        n = len(seq)
        if n <= 1:
            return seq
        mid   = int(n / 2)
        left  = permuteSequence(seq[:mid])
        right = permuteSequence(seq[mid+1:])
        ret   = [seq[mid]]
        while left or right:
            if left:
                ret.append(left.pop(0))
            if right:
                ret.append(right.pop(0))
        return ret

    logging.info('Calculating grid coordinates of SVM parameter')
    cSeq = permuteSequence(rangeF(SVM_CSTART, SVM_CEND, SVM_CSTEP))
    gSeq = permuteSequence(rangeF(SVM_GSTART, SVM_GEND, SVM_GSTEP))

    nC   = float(len(cSeq))
    nG   = float(len(gSeq))
    i    = 0
    j    = 0
    jobs = []
    while i < nC or j < nG:
        if i / nC < j / nG:
            # increase C resolution
            line = []
            for k in range(0, j):
                line.append((cSeq[i], gSeq[k]))
            i = i + 1
            jobs.append(line)
        else:
            # increase g resolution
            line = []
            for k in range(0, i):
                line.append((cSeq[k], gSeq[j]))
            j = j + 1
            jobs.append(line)
    return jobs

class SVMGridWorkerStopToken:
    """Notify the worker to stop"""
    pass

class SVMGridWorker(threading.Thread):
    """Worker thread that calls successive svm-train commands"""

    def __init__(self, name, jobQueue, resultQueue, dataPath):
        threading.Thread.__init__(self)
        self.name        = name
        self.jobQueue    = jobQueue
        self.resultQueue = resultQueue
        self.dataPath    = dataPath

    def run(self):
        """Call svm-train jobs until all the grid coordinates have been explored"""

        while True:
            (c, g) = self.jobQueue.get()
            if c is SVMGridWorkerStopToken:
                self.jobQueue.put((c, g))
                break
            try:
                rate = self.runeOne(2.0 ** c, 2.0 ** g)
                if rate is None:
                    raise RuntimeError(RuntimeError("Got no rate"))
            except:
                # We failed, let others do that and we just quit
                excInfo = sys.exc_info()
                msg     = traceback.format_exception(excInfo[0], excInfo[1], excInfo[2])
                msg     = ''.join(msg)
                logging.warning('[WARNING] Worker %s failed:' % self.name)
                logging.warning(msg)
                self.jobQueue.put((c, g))
                break
            else:
                self.resultQueue.put((self.name, c, g, rate))

    def runeOne(self, c, g):
        """Call a single svm-train job"""

        svmTrain = checkExecutable('svm-train')
        cmd      = [svmTrain,
                    '-c',
                    str(c),
                    '-g',
                    str(g),
                    '-v',
                    str(SVM_FOLD),
                    self.dataPath]
        cmd      = ' '.join(cmd)
        proc     = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        result   = proc.stdout.readlines()
        for line in result:
            if line.find(bytes("Cross",'utf8')) != -1:
                return float(line.split()[-1][:-1])

def doSVMGrid(options, dataPath):
    """Search a parameter grid to optimise the SVM parameters.  This is roughly
    equivalent to the svm-grid script from libsvm."""

    logging.info('Optimising SVM parameters')
    jobs        = calculateSVMGridJobs()
    jobQueue    = Queue.Queue(0)
    resultQueue = Queue.Queue(0)
    for line in jobs:
        for (c, g) in line:
            jobQueue.put((c, g))
    jobQueue._put = jobQueue.queue.appendleft

    for i in range(options.args.nbThreads):
        worker = SVMGridWorker('Worker%03d' % i, jobQueue, resultQueue, dataPath)
        worker.start()

    doneJobs   = {}
    svmOutPath = options.getSVMOutPath()
    resultFile = open(svmOutPath, "w")
    bestRate   = -1
    bestC1     = 0
    bestG1     = 0
    bestC      = 1
    bestG      = 1

    for line in jobs:
        for (c, g) in line:
            while (c, g) not in doneJobs:
                (workerName, c1, g1, rate) = resultQueue.get()
                doneJobs[(c1, g1)]         = rate
                resultFile.write('%f %f %f\n' % (c1, g1, rate))
                if (rate > bestRate) or (rate == bestRate and g1 == bestG1 and c1 < bestC1):
                    bestRate = rate
                    bestC1   = c1
                    bestG1   = g1
                    bestC    = 2.0 ** c1
                    bestG    = 2.0 ** g1
    jobQueue.put((SVMGridWorkerStopToken, None))
    resultFile.close()
    logging.info('Optimal SVM parameters: c=%f, g=%f, rate=%f' % (bestC, bestG, bestRate))
    return bestC, bestG, bestRate

def loadSVMPredictions(path):
    """Loads SVM classifications"""

    handle      = open(path)
    content     = handle.read()
    predictions = content.strip().split('\n')
    handle.close()
    return predictions


def loadBlastClassification(options):
    """Loads previously computed blast classifications"""

    blastSort      = options.getBlastSortPath()
    handle         = open(blastSort)
    classification = {}
    n              = 0
    for line in handle:
        line = line.strip()
        if not line:
            continue
        fields                 = line.split()
        seqId                  = fields[0]
        seqCode                = fields[1]
        classification[seqId]  = seqCode
        n                     += 1
    logging.info('Parsed %d blast classifications' % n)
    return classification

def writeOutput(options, predictions, blastClassification, fastaPath, fastaName, prefix1, prefix2):
    """Write the final results"""

    logging.info('Writing final output files')
    # Create output directory and open output files
    size          = len(predictions)
    hostResults   = prefix1 + '_' + fastaName
    symbResults   = prefix2 + '_' + fastaName
    if options.args.outDir:
        outFolder = os.path.abspath(options.args.outDir)
        if not os.path.isdir(outFolder):
            os.makedirs(outFolder)
        hostHandle  = open(os.path.join(outFolder, hostResults), "w")
        symbHandle  = open(os.path.join(outFolder, symbResults), "w")
    else:
        hostHandle  = open(hostResults, "w")
        symbHandle  = open(symbResults, "w")

    # Iterate over sequences and save them according to their classification
    tot           = 0
    seqIdx        = 0
    disagreements = 0
    for name, seq in iterFasta(fastaPath):
        tot      += 1
        name      = name.split()[0]
        blastCode = blastClassification.get(name, UNKNOWN_CODE)
        if predictions[seqIdx] == blastCode or blastCode == UNKNOWN_CODE:
            if predictions[seqIdx] == HOST_CODE_STR:
                hostHandle.write('>%s\n%s\n' % (name, seq))
            elif predictions[seqIdx] == SYMB_CODE_STR:
                symbHandle.write('>%s\n%s\n' % (name, seq))
        # Keep blast classification if it disagrees with SVM
        elif predictions[seqIdx] != blastCode:
            disagreements += 1
            if blastCode == HOST_CODE_STR:
                hostHandle.write('>%s\n%s\n' % (name, seq))
            elif blastCode == SYMB_CODE_STR:
                symbHandle.write('>%s\n%s\n' % (name, seq))
        seqIdx += 1
        if seqIdx > size:
            logging.warning('[WARNING] Found more sequences than prediction. This may be caused by dupplicated sequence names.')
            break

    logging.info('Classified %d sequences, found %d disagreements between SVM and blast.' % (tot, disagreements))
    hostHandle.close()
    symbHandle.close()

def predictSVM(options, blastClassification, kmerTrain, kmerTest):
    """Final SVM predictions and combination with the blast results"""

    logging.info('Predicting with SVM optimal parameters')
    svmPredict = checkExecutable('svm-predict')
    svmScale   = checkExecutable('svm-scale')
    kmerTrain  = os.path.join(options.tempDir, kmerTrain)
    modelFile  = kmerTrain + '.model'
    rangeFile  = kmerTrain + '.range'
    resultLog  = kmerTrain + '_svmPredict.log'
    fastaPath  = options.args.queries
    fastaName  = os.path.basename(fastaPath)
    kmerScale  = os.path.join(options.tempDir, fastaName + '.scaled')
    kmerPred   = os.path.join(options.tempDir, fastaName + '.pred')
    kmerFile   = os.path.join(options.tempDir, fastaName + '.kmers')
    computeKmers(options, options.args.queries, fastaName + '.kmers', HOST_CODE, "w")
    #SVM_Scale
    scaleCmd   = [svmScale,
                  '-r',
                  rangeFile,
                  kmerFile,
                  '>',
                  kmerScale]
    scaleCmd   = ' '.join(scaleCmd)
    retCode    = subprocess.call(scaleCmd, shell=True)
    if not retCode == 0:
        logging.error('[ERROR] Please check inputs. svm-scale not executed or exit with error.')
        sys.exit(1)
    #SVM_predict
    predictCmd = [svmPredict,
                  kmerScale,
                  modelFile,
                  kmerPred]
    predictCmd   = ' '.join(predictCmd)
    resultHandle = open(resultLog, "w")
    retCode      = subprocess.call(predictCmd, shell=True, stdout=resultHandle)
    resultHandle.close()
    if not retCode == 0:
        logging.error('[ERROR] Please check inputs. svm-predict not executed or exit with error.')
        sys.exit(1)
    #parse_Prediction
    predictions = loadSVMPredictions(kmerPred)
    writeOutput(options, predictions, blastClassification, fastaPath, fastaName, HOST_NAME, SYMB_NAME)

def checkExecutable(program):
    """Check whether a program is installed and executable"""

    # First check in $PATH
    path = os.getenv('PATH')
    for d in path.split(os.path.pathsep):
        exe = os.path.join(d, program)
        if os.path.exists(exe) and os.access(exe, os.X_OK):
            return exe
    # Then check in the subdirectory
    root = os.path.dirname(os.path.abspath(sys.argv[0]))
    exe  = os.path.join(root, BINARIES_DIR, program)
    if os.path.exists(exe) and os.access(exe, os.X_OK):
        return exe

def mainArgs():
    """Process command-line arguments"""

    parser = argparse.ArgumentParser(description='Perform SVM Classification of Host and Symbiont (or Parasite) Sequences')
    parser.add_argument('queries',
                        help='The input queries sequences')

    parser.add_argument('-A',
                        '--species1',
                        type=str,
                        help='Reference sequences for the first species')
    parser.add_argument('-B',
                        '--species2',
                        type=str,
                        help='Reference sequences for the second species')
    parser.add_argument('-b',
                        '--blastResults',
                        type=str,
                        help='Blast results obtained')
    parser.add_argument('-T',
                        '--blastType',
                        default='blastx',
                        choices=('blastx', 'blastn', 'tblastx'),
                        help='Type of blast search to be performed')
    parser.add_argument('-p',
                        '--nbThreads',
                        type=int,
                        default='1',
                        help='Number of threads to run the BLAST search and SVM')
    parser.add_argument('-e',
                        '--maxBestEvalue',
                        type=float,
                        default='1e-20',
                        help='Maximum e-value')
    parser.add_argument('-n',
                        '--numberOfSeq',
                        type=int,
                        default='1000',
                        help='Maximum number of sequences for training and testing')
    parser.add_argument('-s',
                        '--minSeqSize',
                        type=int,
                        default='0',
                        help='Minimum sequence size for training and testing')
    parser.add_argument('-c',
                        '--minWordSize',
                        type=int,
                        default='1',
                        help='Minimum value of DNA word length')
    parser.add_argument('-k',
                        '--maxWordSize',
                        type=int,
                        default='4',
                        help='Maxmimum value of DNA word length')
    parser.add_argument('-r',
                        '--bothStrands',
                        action='store_true',
                        help='Compute kmers for the forward and reverse strands')
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='Turn Verbose mode on?')
    parser.add_argument('-t',
                        '--tempDir',
                        type=str,
                        default='temp',
                        help='Location (prefix) of the temporary directory')
    parser.add_argument('-o',
                        '--outDir',
                        default='',
                        help='Name of optional output directory')
    parser.add_argument('-X',
                        '--clearTemp',
                        action='store_true',
                        help='Clear all temporary data upon completion')
    parser.add_argument('-z',
                        '--stopAfter',
                        choices=('db', 'runBlast', 'parseBlast', 'kmers', 'SVM'),
                        help='Optional exit upon completion of stage.')
    parser.add_argument('-R',
                        '--restart',
                        action='store_true',
                        help='Continue process from last exit stage')
    parser.add_argument('-V',
                        '--version',
                        action='version',
                        version='%(prog)s ' + PSYTRANS_VERSION)
    args = parser.parse_args()
    if args.minWordSize > args.maxWordSize:
        logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
        sys.exit(1)
    return args

def main():
    """Banzai !!!"""

    # Initialise RNG
    random.seed(1234)

    # Get the options
    args    = mainArgs()
    options = PsyTransOptions(args)

    # Setup loging to a file and to the console
    logFormat = "%(asctime)s - %(funcName)s - %(message)s"
    logging.basicConfig(level=logging.INFO,
                        format=logFormat,
                        filename=options.getLogPath(), filemode="w")
    console   = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter(logFormat)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    logging.info("Arguments parsed. Starting...")

    programList = [args.blastType, 'svm-train', 'svm-scale', 'svm-predict']
    for program in programList:
        if not checkExecutable(program):
            logging.warning('[WARNING] %s program not found !!! This may cause problems later on.' % program)

    # Check that we are either starting from blast results or references sequences
    if not args.blastResults and not (args.species1 and args.species2):
        logging.error('[ERROR] Either provide references for both species or blast results')
        sys.exit(1)

    # Start from the input sequences
    if not args.blastResults:
        #Step 1: make the database
        if not (args.restart and options.checkPoint("writeDatabase.done")):
            writeDatabase(options, options.getFastaDbPath())
            if checkExecutable('makeblastdb'):
                makeDB(options)
            else:
                logging.error('[ERROR] makeblastdb not found. Exiting')
                sys.exit(1)
        if args.stopAfter == 'db':
            logging.info('Stop after "db" requested, exiting now')
            sys.exit(0)
        #Step 2: run the blast searches
        if not (args.restart and options.checkPoint("runBlast.done")):
            if checkExecutable('blastx'):
                runBlastThreads(options)
            else:
                logging.error('[ERROR] blastx not found. Exiting')
                sys.exit(0)
        if args.stopAfter == 'runBlast':
            logging.info('Stop after "runBlast" requested, exiting now')
            sys.exit(0)
    # Start from the user-provided blast results
    elif not os.path.exists(args.blastResults):
        logging.error('[ERROR] Could not find user-provided blast results (%s). Exiting' % args.blastResults)
        sys.exit(1)

    #Step 3: process the blast results
    if not (args.restart and options.checkPoint("parseBlast.done")):
        processBlastResults(options)
    if args.stopAfter == 'parseBlast':
        logging.info('Stop after "parseBlast" requested, exiting now')
        sys.exit(0)

    #Step 4: Kmer preparation
    kmerTrainPath = options.getTrainPath()
    kmerTestPath  = options.getTestPath()
    if not (args.restart and options.checkPoint("kmers.done")):
        prepareTrainingKmers(options, kmerTrainPath, kmerTestPath)
    if args.stopAfter == 'kmers':
        logging.info('Stop after "kmers" requested, exiting now')
        sys.exit(0)

    #Step 5: SVM training and testing
    if not (args.restart and options.checkPoint("svm.done")):
        if checkExecutable('svm-train') and checkExecutable('svm-scale') and checkExecutable('svm-predict'):
            doSVMEasy(options, kmerTrainPath, kmerTestPath)
        else:
            logging.error('[ERROR] Failed to find some of the libsvm commands. Make sure that svm-train, svm-scale and svm-predict are installed.')
            sys.exit(1)
    if args.stopAfter == 'SVM':
        logging.info('Stop after "SVM" requested, exiting now')
        sys.exit(0)

    #Step 6: final classification
    blastClassification = loadBlastClassification(options)
    predictSVM(options, blastClassification, kmerTrainPath, kmerTestPath)
    logging.info("SVM classification completed successfully.")

    if args.clearTemp:
        shutil.rmtree(options.tempDir)

if __name__ == '__main__':
    main()

# vim:ts=4:sw=4:sts=4:et:ai:
