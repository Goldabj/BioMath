'''
Created on Mar 21, 2017

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
from seaborn.utils import sig_stars


def SequenceExamples(): 
    dna = Seq("TATGACCAGGGAGCTCTGGAT", IUPAC.unambiguous_dna)
    reverse_complement = dna.reverse_complement()
    
    for i in range (0, 3) : 
        readingFrame = dna[i:]
        protein = readingFrame.translate(stop_symbol = "*")
        print("protein " + str(i) + " = "+ str(protein))
        
    for i in range (0, 3) : 
        readingFrame = reverse_complement[i:]
        protein = readingFrame.translate(stop_symbol = "*")
        print("protein " + str(i+3) + " = " + str(protein))

def FindSequence() :
    dna = Seq("AACTGATGGTACCCTACAGGAATTGTACCCGATTTTGAGGTC")
    print(dna.transcribe())
    
    print('\n')
    
    index = dna.find("ATG")
    print(dna[index:].translate())
    
    
def readBlosum62AlignmentScores() : 
    blosum62 = pd.read_csv("file:///C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/insulin_101_AA_counts.csv")
    blosum62.index = blosum62.columns
    

    seq1 = "LALWGP"
    seq2 = "IWSVGT"
    
    for s1 in seq2:
        for s2 in seq1:
            print(str(blosum62.get_value(s1, s2)) + " ", end=' ')
        print('\n')
    
def probablitiyMutationExample():
    # subsMatrix = round( 10 * log(P(S n R)/(P(S) * P(R)))) for all S & R 
    
    counts_raw = np.array([[5, 1, 0], [0, 3, 3], [0, 1, 3]])
    counts_raw = pd.DataFrame(counts_raw, index=["S", "N", "R"], columns=["S", "N", "R"])
    counts = counts_raw + 1
    
    # intersection probablitiy matrix 
    intersection_probs = counts /(counts.sum().sum())
    
    # denomonator problility matrix
    # counts.sum() = individual counts , np.outer() joins the two matrixes on all combonations of each 
    denom = np.outer(counts.sum()/counts.sum().sum(), counts.sum()/counts.sum().sum())
    
    #subs matrix - round and multimply by 10
    subsMatrix = round(10 * np.log10(intersection_probs / denom))
    
    print(subsMatrix)
    
def computeMutationMatrix() :
    countsMatrix = np.array([[5, 1, 0], [0, 3, 3], [0, 1, 3]])
    
    countsMatrix = pd.DataFrame(countsMatrix, index=["S", "N", "R"], columns=["S", "N", "R"]) # give the array index and columsn
    
    countsMatrix = countsMatrix + 1 # pseudo counts
    
    countsMatrix = 0.5 * (countsMatrix + countsMatrix.T)  # makes the counts symmetric 
    
    probsMatrix = countsMatrix/countsMatrix.sum().sum()  # P(A n B)/P(A)
    
    p = probsMatrix.sum() # probs of S, N, and R 
    
    MutationMatrix = probsMatrix.div(p, axis=0) # divide each element by p
    
    1-np.diag(MutationMatrix) # the mutation rates. Diagolnals are the probablitity that they do not mutate
    
    mutability = pd.Series(1-np.diag(MutationMatrix), index=MutationMatrix.columns) #mutation rates in a series 
    
    totalMutationRate = np.dot(mutability, p) # mutationMatrix  * p(start) 
    
    # ------------- To scale to 1% mutation matrix -------------------------------------------- #
    Mii = np.diag(mutability) # creates a diagonal matrix of the mutabilities 
    Mij = MutationMatrix - Mii # removes the diagonals from mutation matrix
    
    # get the 1% mutation matrix: rate = p(mutation | letter) * p(letter) 
    # multiple = 1/rate -> .01 = z * (p(mut | letter) * p(letter)) 
    # z = 0.01/(np.dot(Mij , p).sum())
    # Mij1 = z* Mij ---- is the mutation matrix adjusted to 1% -- however the rows do not add to 1
    # then caluculate what the diagonals should be to get MutationMatrix%1 Mij + np.diag(1-Mij1.sum(axis=1))
    
    
def alignmentStats5_2():
    # generate random sequences
    '''seq1 = ''.join(np.random.choice(["N", "S", "R"], p=[0.1, 0.2, 0.7], size=10)) # letters, distributions, and size
    
    # compute scoring matrix 
    aa_letters = list("A...") # the amino letters
    aa_counts = pd.read_csv("file...") # aacounts index_col = 0
    
    aa_probs = aa_counts/aa_counts.sum()
    
    # calculate random scores a bunch of times
    seq1 = randomSeq(aa_letters, aa_probs.sequeze(), 100)
    seq2 = randomSeq(aa_letters, aa_probs.sequeze(), 100)
    score = computeScore(seq1, seq2, blosum62)
    # plot the distribution of scores of random sequences 
    scr = pd.Sieres(allScores)
    scr.plot.density()'''
    
    
def pscorestuff():
    b62 = matlist.blosum62
    b62 = b62.drop(['B', 'Z']) # there are more columns to drop
    b62 = b62.drop(['B', 'Z']) # B, Z , X , and *a
    # drop the other row
    b62.mean().mean() # without proper distribution 
    
    aa_counts = pd.read_csv("file:///C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/insulin_101_AA_counts.csv", index_col=0)
    
    aa_probs = aa_counts/aa_counts.sum()
    
    aa_probs_matchup = np.outer(aa_probs, aa_probs)
    
    mean = (b62 * aa_probs_matchup).sum()
    
    meanS = 110* mean
    
    sig = np.sqrt((b62 -mean) * (b62*mean) *np.outer(aa_probs, aa_probs)).sum().sum()
   
    sigS = np.sqrt(110) * sig
    
    z_score = (25-meanS) / sigS
    
    
    
    
    
