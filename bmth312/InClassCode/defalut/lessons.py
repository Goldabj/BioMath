'''
Created on Apr 6, 2017

@author: Goldacbj
'''
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import scipy.stats as stats

def computePam1() :
    countMatrix = pd.read_csv("file:///C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/PAMMartrices_countMatrix85.csv")
    
    totalCount = countMatrix.sum().sum()
    
    probsMatrix = countMatrix / totalCount
    
    mutationMatrix = probsMatrix.div(probsMatrix.sum(), axis=1)
    
    mutability = np.diag(mutationMatrix)
    
    
def lesson417() :
    human_insulin = getSequenceFromFile("sp|P01308|INS_HUMAN")
    zScore = computeZscoreForAlignment(value=25, referenceSeq=human_insulin, lenOfRandomSequence=110, globalAlignment=True)
    pScore = computePscoreFromZscore(zScore)
    print("1-PartB: z-score for human to snake alignment: " + str(zScore))
    print("         p-score for human to snake alignment: " + str(pScore))
    print("")
    
    

# ----------- helpers ---------------------------------------
    
def generateRandomSequences(letters=None, distribution=None, sequenceSize=10, numOfSequences=100):
    ''' returns a list of random sequences of size specified, distributions should be a map from letter to its 
    probablility'''
    dist = list()
    if (letters == None) :
        letters = list("ARNDCQEGHILKMFPSTWYV")
    if (distribution == None) : 
        dist = list(1 / (len(letters)) for i in range(0, len(letters)))
    else :
        # make sure the distributions are ordered correctly
        for letter in letters: 
            dist.append(distribution[letter])
    
    return list(''.join(np.random.choice(a=letters, p=dist, size=sequenceSize)) for i in range(0, numOfSequences))

def getInsulinAlignmentScore(proteinKey1, proteinKey2, globalAlign=True):
    ''' gets the insulin protein keys from file, then aligns them and returns the best score'''
    subsMatrix = matlist.blosum62
    seq1 = getSequenceFromFile(proteinKey1)
    seq2 = getSequenceFromFile(proteinKey2)
    
    score = 0
    if (globalAlign) :
        score = pairwise2.align.globalds(seq1, seq2, subsMatrix, -11, -1, score_only=True)
    else : 
        score = pairwise2.align.localds(seq1, seq2, subsMatrix, -11, -1, score_only=True)
    
    return score
    
def getSequenceFromFile(proteinKey):
    ''' returns the sequence from the insulin file'''    
    other_insulins = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/insulin_101_uniprot.fasta.txt", "fasta"))
    for record in other_insulins : 
        if (record.id == proteinKey) :
            return record.seq
    
    return None

def computeZscoreForAlignment(value, referenceSeq, lenOfRandomSequence, globalAlignment=True):
    ''' computes the z score for the alignment value by aligning the referenceSequence to 1000 random sequences'''
    randomScores = computeRandomAlignmentScores(referenceSeq, lenOfRandomSequence, globalAlignment)
    print("random alignment scores mean: "+ str(np.mean(randomScores)))
    print("random alignment scores std: " + str(np.std(randomScores)))
    z_score = (value - np.mean(randomScores))/np.std(randomScores)
    
    return z_score; 
    
def computePscoreFromZscore(z_score) : 
    ''' computes the p score for the alignment value by aligning the referenceSequence to 1000 random sequences'''
    p_score = 1 - stats.norm.cdf(z_score);
    return p_score;
    
def computeRandomAlignmentScores(seq1, lenOfRandomSequence, globalAlignment=True):
    ''' computes random alignment scores and returns a list of scores'''
    subsMatrix = matlist.blosum62
    randomScores = list()
    distribution = computeLetterDistribution(seq1);
    
    randomSequences = generateRandomSequences(distribution=distribution, sequenceSize=lenOfRandomSequence, numOfSequences=1000)
    
    for seq in randomSequences : 
        score = 0
        if (globalAlignment) : 
            score = pairwise2.align.globalds(seq1, seq, subsMatrix, -11, -1, score_only=True)
        else : 
            score = pairwise2.align.localds(seq1, seq, subsMatrix, -11, -1, score_only=True)
        randomScores.append(score)
        
    return randomScores

def computeLetterDistribution(seq):
    ''' computes the letter distribution probabilites from the sequence and returns a map of these probs'''
    counts = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0}
    for letter in seq: 
        prev_count = counts[letter]
        counts[letter] = prev_count + 1
    
    sum = 0
    for key in counts.keys(): 
        sum = sum + counts[key]
    
    for key in counts.keys(): 
        counts[key] = counts[key]/sum
        
    return counts

lesson417()