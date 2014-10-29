#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

###############################################################################
###############################################################################
##                                                                           ##
##  MIS-203                                                                  ##
##  Clustering Lecture   (April 12, 2005)                                    ##
##  Programming Assignment                                                   ##
##                                                                           ##
##  Author: Scott Pegg                                                       ##
##                                                                           ##
##                                                                           ##
##  The abbreviated instructions:                                            ##
##    You will be given a set of enzyme active sites in PDB format           ##
##    (1) Implement a similarity metric                                      ##
##    (2) Implement a clustering method based on a partitioning algorithm    ##
##    (3) Implement a clustering method based on a hierarchical algorithm    ##
##    (4) Answer the questions given in the homework assignment              ##
##                                                                           ##
##  Please read the full instructions from the course website _before_ you   ##
##  start this assignment!                                                   ##
##                                                                           ##
###############################################################################
###############################################################################



from string import *
from math import *
import numpy as np
import sys, os
import copy
import glob
import pdb
import random
from collections import Counter
from Bio import SeqIO

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) != 7:
  print "Usage: smith-waterman.py <substution matrix> <compfile> <output file> <gap penalty> <extension penalty>"
  sys.exit(0)


###############################################################################
#                                                                             #
# A simple class for a sequence                                               #

class Sequence:
    def read_fasta(self, file):
        f = open(file)
        c = SeqIO.parse(f,'fasta')
        i = 0
        for b in c:
            i = b
        #pdb.set_trace()
        return (i.id, str(i.seq))
           
        
    def __init__(self, fileLocation):
        self.fileLocation = fileLocation
        (self.name, self.sequence) = self.read_fasta(fileLocation)
    
    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return str(self.name) + str(self.sequence)

###############################################################################



###############################################################################
#                                                                             #
# A simple class for an amino acid residue                                    #

class AlignmentMatrix:
    def calculateGapScore(self, pointers, row, col, gap_pen, ext_pen, movingTo):
        if pointers[row][col] == 2 and movingTo == 'down':
            #print 'YES!!!'
            return ext_pen
        if pointers[row][col] == 1 and movingTo == 'right':
            #print 'YESSSSSS'
            return ext_pen
        else:
            return gap_pen
    def calculateAlignments(self, seq1, seq2, subsMatrix, gap_pen, ext_pen):
        # Simplify gap penalties

        best_score = 0
        best_i = 0
        best_j = 0
        
        # Add in a gap to start
        if seq1.sequence[0] != '*':
            seq1.sequence = '*' + seq1.sequence

        if seq2.sequence[0] != '*':
            seq2.sequence = '*' + seq2.sequence
            
        alignmentScores = [ [ 0 for j in range(len(seq2.sequence)) ] for i in range(len(seq1.sequence)) ]
        pointer = [ [ 0 for j in range(len(seq2.sequence)) ] for i in range(len(seq1.sequence)) ]
        # pdb.set_trace()
        # Down (rows)
        for i in range(1,len(seq1.sequence)):
            # Right (cols)
            for j in range(1,len(seq2.sequence)):
                # print i, j
                scoreRight = alignmentScores[i][j-1] - self.calculateGapScore(pointer, i, j-1, gap_pen, ext_pen, 'right')
                scoreDown = alignmentScores[i-1][j] - self.calculateGapScore(pointer, i-1, j, gap_pen, ext_pen, 'down')
                scoreDiagonal = alignmentScores[i-1][j-1] + subsMatrix.getScore(seq1.sequence[i],seq2.sequence[j])
                possibleScores = [scoreRight, scoreDown, scoreDiagonal, 0]
                alignmentScores[i][j] = max(possibleScores)
                if alignmentScores[i][j]==0:
                  pointer[i][j]=0; #0 means end of the path
                if alignmentScores[i][j]==scoreRight:
                  pointer[i][j]=1; #1 means trace came from the left
                if alignmentScores[i][j]==scoreDown:
                  pointer[i][j]=2; #2 means trace came from above
                if alignmentScores[i][j]==scoreDiagonal:
                  pointer[i][j]=3; #3 means trace came from up and left diagonal
                if alignmentScores[i][j] > best_score:
                    best_score = alignmentScores[i][j]
                    best_i = i
                    best_j = j
                    
        #best_score = float(best_score)/float(min(len(seq1.sequence),len(seq2.sequence)))
        #print alignmentScores, pointer, best_score, best_i, best_j       
        return alignmentScores, pointer, best_score, best_i, best_j
  
    def __init__(self, seq1, seq2, subsMatrix, gap_pen, ext_pen):
        self.alignments, self.pointer, self.best_score, self.best_i, self.best_j = self.calculateAlignments(seq1, seq2, subsMatrix, gap_pen, ext_pen)
  
    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.alignments

###############################################################################



###############################################################################
#                                                                             #
# A simple class for an active site                                           #

class SubstitutionMatrix:
    
    def getScore(self, l1, l2):
        l1_index = self.index[l1]
        l2_index = self.index[l2]
        return self.scores[l1_index][l2_index]
        
    def randomize(self, r):
        i = [random.randrange(0,len(self.index)) for k in range(15)]
        j = [random.randrange(0,len(self.index)) for l in range(15)]
        for m in range(15):
            move = random.gauss(0,5)
            self.scores[i[m]][j[m]] += move
            self.scores[j[m]][i[m]] += move
        
    def makenice(self):
        # score = 0.0
#         for a in range(len(self.index)):
#             for b in range(len(self.index)):
#                 score += self.scores[a][b]
#         return score
        outstring = ''
        for i in range(len(self.index)):
            outstring += self.index.keys()[i] + ' '
        outstring += '\n'
        for i in range(len(self.index)):
            for j in range(len(self.index)):
                outstring += str(self.scores[i][j]) + ' '
            outstring += '\n'
        return outstring
        
    def buildMatrix(self, file):
        lines = open(file).readlines()
        FirstTrueLinePassed = False
        index = {}
        # Hack, since we dont know how many letters will show up on the first
        substutionScores = [[]]
        row = -1
        for line in lines:
            print line
            print len(line.split())
            if line[0] != '#':
                if FirstTrueLinePassed:
                    # Build the matrix
                    row += 1
                    substutionScores[row] = [float(num) for num in line.split()]
                else:
                    # Set up the index on the first pass
                    FirstTrueLinePassed = True
                    letters = line.split()
                    print letters
                    for i, l in enumerate(letters):
                        index[l] = i
                    size = len(letters)
                    substutionScores = [ [ 0 for i in range(size) ] for j in range(size) ]
        return index, substutionScores
                
  
    def __init__(self, name):
        if isinstance(name, str):
            fileLocation = name
            self.fileLocation = fileLocation
            self.index, self.scores = self.buildMatrix(fileLocation)
        elif isinstance(name, SubstitutionMatrix):
            anotherself = name
            self.fileLocation = copy.deepcopy(anotherself.fileLocation)
            self.index = copy.deepcopy(anotherself.index)
            self.scores = copy.deepcopy(anotherself.scores)
        else:
            print 'FAIL'
  
        # Overload the __repr__ operator to make printing simpler.
        def __repr__(self):
            outstring = ''
            for i in range(len(self.index)):
                outstring += self.index.keys()[i] + ' '
            outstring += '\n'
            for i in range(len(self.index)):
                for j in range(len(self.index)):
                    outstring += str(self.scores[i][j]) + ' '
                outstring += '\n'
            return outstring
            
        def __str__(self):
            outstring = ''
            for i in range(len(self.index)):
                outstring += self.index.keys()[i] + ' '
            outstring += '\n'
            for i in range(len(self.index)):
                for j in range(len(self.index)):
                    outstring += str(self.scores[i][j]) + ' '
                outstring += '\n'
            return outstring

###############################################################################



###############################################################################
#                                                                             #
# Read in all of the active sites from the given directory.                   #
#                                                                             #
# Input: directory                                                            #
# Output: list of ActiveSite instances                                        #

def read_sequences(dir):
  
  files = glob.glob(dir + '/*.fa')
  
  sequences = []
  for file in files:
    sequence = Sequence(file)
    
    sequences.append(sequence)
  return sequences
  
def read_comparisons(compFile):
    f = open(compFile)
    comparisons = []
    lines = f.readlines()
    for line in lines:
        files = line.split()
        seq1 = Sequence(files[0])
        seq2 = Sequence(files[1])
        comp = (seq1, seq2)
        comparisons.append(comp)
    return comparisons

###############################################################################
#                                                                             #
# Top Level                                                                   #
# "Usage: smith-waterman.py <substution matrix> <comparison list> <output file> <gap penalty> <extension penalty>"

def trackBack(pointers, seq1, seq2, i, j):
    '''Tracks back to create the aligned sequence pair'''
    alignedSeq1 = ''
    alignedSeq2 = ''
 
    #print 'Trackback will start at %d,%d' % (i, j)
    while pointers[i][j] != 0:
        #print 'Trackback: Working on %d,%d = %d' % (i, j, pointers[i][j])
        #if i == 0 or j == 0:
            #print 'FUCK: 0 hate non-terminating condition'
            #break
        if pointers[i][j] == 1:
            # print 'Trackback: Delete'
            alignedSeq1 = seq1[i - 1] + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif pointers[i][j] == 2:
            # print 'Trackback: Insert'
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = seq2[j - 1] + alignedSeq2
            j = j - 1
        elif pointers[i][j] == 3:
            # print 'Trackback: Match'
            alignedSeq1 = seq1[i - 1] + alignedSeq1
            alignedSeq2 = seq2[j - 1] + alignedSeq2
            i = i - 1
            j = j - 1
        # print 'Trackback: Done %d,%d' % (i, j)
 
    return (alignedSeq1, alignedSeq2, i, j)

###############################################################################

###############################################################################
#                                                                             #
# Top Level                                                                   #
# "Usage: smith-waterman.py <substution matrix> <comparison list> <output file> <gap penalty> <extension penalty>"

def getROC(pointers, seq1, seq2, i, j):
    '''Tracks back to create the aligned sequence pair'''
    alignedSeq1 = ''
    alignedSeq2 = ''
 
    #print 'Trackback will start at %d,%d' % (i, j)
    while pointers[i][j] != 0:
        #print 'Trackback: Working on %d,%d = %d' % (i, j, pointers[i][j])
        #if i == 0 or j == 0:
            #print 'FUCK: 0 hate non-terminating condition'
            #break
        if pointers[i][j] == 1:
            # print 'Trackback: Delete'
            alignedSeq1 = seq1[i - 1] + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i = i - 1
        elif pointers[i][j] == 2:
            # print 'Trackback: Insert'
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = seq2[j - 1] + alignedSeq2
            j = j - 1
        elif pointers[i][j] == 3:
            # print 'Trackback: Match'
            alignedSeq1 = seq1[i - 1] + alignedSeq1
            alignedSeq2 = seq2[j - 1] + alignedSeq2
            i = i - 1
            j = j - 1
        # print 'Trackback: Done %d,%d' % (i, j)
 
    return (alignedSeq1, alignedSeq2, i, j)

###############################################################################

# Compare items in a list
if sys.argv[1][0:2] == '-C':
    comparisons = read_comparisons(sys.argv[3])
    subsMatrix = SubstitutionMatrix(sys.argv[2])
    outputFile = sys.argv[4]
    gap_pen = int(sys.argv[5])
    ext_pen = int(sys.argv[6])
    fout = open(outputFile, 'w')
    for comp in comparisons:
        seq1 = comp[0]
        seq2 = comp[1]
        # pdb.set_trace()
        a = AlignmentMatrix(seq1, seq2, subsMatrix, gap_pen, ext_pen)
        s1, s2, i, j = trackBack(a.pointer, seq1.sequence, seq2.sequence, a.best_i, a.best_j)
        print a.best_score
        print seq1.fileLocation
        print seq2.fileLocation
        print s1
        print s2
        outString = str(seq1.fileLocation) + '\t' + str(seq2.fileLocation) + '\t' + str(a.best_score) + '\t' +str(sys.argv[2])
        print>>fout, outString
            

# Optimize a matrix
if sys.argv[1][0:2] == '-O':
    comparisonFiles = sys.argv[3].split("::")
    posComparisons = read_comparisons(comparisonFiles[0])
    negComparisons = read_comparisons(comparisonFiles[1])
    outputFiles = sys.argv[4].split("::")
    posOutputFile = outputFiles[0]
    negOutputFile = outputFiles[1]
    gap_pen = int(sys.argv[5])
    ext_pen = int(sys.argv[6])
    converged = False
    bestScore = 0
    bestSubsMatrix = SubstitutionMatrix(sys.argv[2])
    currentSubsMatrix = SubstitutionMatrix(sys.argv[2])
    bestScoreEver = 0
    bestSubsMatrixEver = SubstitutionMatrix(sys.argv[2])
    itter = 0
    
    while not converged:
        #print 'Begin While loop' + currentSubsMatrix.makenice()
        print itter
        itter += 1
        posfout = open(posOutputFile, 'w')
        negfout = open(negOutputFile, 'w')
        for comp in posComparisons:
            seq1 = comp[0]
            seq2 = comp[1]
            a = AlignmentMatrix(seq1, seq2, currentSubsMatrix, gap_pen, ext_pen)
            s1, s2, i, j = trackBack(a.pointer, seq1.sequence, seq2.sequence, a.best_i, a.best_j)
            outString = str(seq1.fileLocation) + '\t' + str(seq2.fileLocation) + '\t' + str(a.best_score) + '\t' +str(sys.argv[2])
            print>>posfout, outString
        
        for comp in negComparisons:
            seq1 = comp[0]
            seq2 = comp[1]
            a = AlignmentMatrix(seq1, seq2, currentSubsMatrix, gap_pen, ext_pen)
            s1, s2, i, j = trackBack(a.pointer, seq1.sequence, seq2.sequence, a.best_i, a.best_j)
            outString = str(seq1.fileLocation) + '\t' + str(seq2.fileLocation) + '\t' + str(a.best_score) + '\t' +str(sys.argv[2])
            print>>negfout, outString
        posfout.close()
        negfout.close()
        
        stripPosCommand = 'cat posOut | cut -f 3 | sort -n -r > pos'
        os.system(stripPosCommand)
        stripNegCommand = 'cat negOut | cut -f 3 | sort -n -r > neg'
        os.system(stripNegCommand)
        getScoreCommand = '../ucsf-roc/a.out pos neg temp'
        
        os.system(getScoreCommand)
        
        lines = open('temp-rocstats', 'r').readlines()
        score = 0.0
        for i in range(4):
            score += float(lines[i].split(' ')[1])
        score = score / 100.0
        print '\n\n\n\nBefore'
        print itter
        print bestScore
        print bestSubsMatrix.makenice()
        if score > bestScore:
            #sys.exit(0)
            bestScore = score
            del bestSubsMatrix
            bestSubsMatrix = SubstitutionMatrix(currentSubsMatrix)
            if score > bestScoreEver:
                bestScoreEver = score
                bestSubsMatrixEver = SubstitutionMatrix(currentSubsMatrix)
            print 'this was better, so i accepted'
        elif score == bestScore:
            print 'i think we are finished'
            # pdb.set_trace()
        elif score < bestScore:
            if random.randrange(1,101) > 90:
                print 'this was worse, but i accepted it anyway'
                bestScore = score
                del bestSubsMatrix
                bestSubsMatrix = SubstitutionMatrix(currentSubsMatrix)
        
        print '\n\n\n\nAfter'
        print itter
        print bestScore
        print bestSubsMatrix.makenice()
        print '\n\n\n\n'
        
        # Randomize the subsmatrix now'
        
        currentSubsMatrix = SubstitutionMatrix(bestSubsMatrix)
        #print 'Before' + currentSubsMatrix.makenice()
        currentSubsMatrix.randomize(2)
        #print 'After' + currentSubsMatrix.makenice()
        bfile = open('bestMatrix', 'w')
        outstring = bestSubsMatrixEver.makenice()
        print outstring
        print bestScoreEver
        print>>bfile, outstring
        bfile.close()
        if itter > 3:
            converged = True
        
        
        
        
        



###############################################################################


