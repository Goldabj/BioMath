'''
Created on Mar 23, 2017

@author: Goldacbj
'''
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np



def interactiveInputAlignment():
    
    seq1 = input("enter protein sequence 1: ")
    seq2 = input("enter protein sequence 2: ")
    
    protein1 = Seq(seq1, IUPAC.extended_protein)
    protein2 = Seq(seq2, IUPAC.extended_protein)
    
    
    blosum62 = pd.read_csv("file:///C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Blosum62.csv")
    blosum62.index = blosum62.columns
    
    _alignment = computeAlignment(protein1, protein2, blosum62, -2)
    
    print(str(_alignment[0]))
    print(str(_alignment[1]) + "\n")
    print("score: ")
    print(str(computeAlignmentScore(alignment=_alignment, gapPentalty=-2, scoringMatrix=blosum62)))
    
    
    

def computeAlignment(protein1, protein2, matchUpMatrix, gapPenalty=-2) :
    
    scoringMatrix = _computeScoringMatrix(matchUpMatrix, protein1, protein2)
    
    print("Scoring Matrix: " + str(scoringMatrix))
    
    optimalPathMatrix = _computeOptimalPathMatrix(scoringMatrix, gapPenalty)
    
    print("optimalPathMatrix: " + str(optimalPathMatrix))
    
    alignment = _backTrackForAlignment(optimalPathMatrix, protein1, protein2)
    
    return alignment

# this method computes the scoring matrix - the dot plot matrix of the 
# two proteins using the scoring matrix passed in
def _computeScoringMatrix(matrix, protein1, protein2):
    scoringMatrix = np.zeros(shape=(len(protein2), len(protein1)))
    
    for i, letter1 in enumerate(protein1) :
        for j, letter2 in enumerate(protein2) : 
            score = matrix.get_value(letter1, letter2)
            scoringMatrix.itemset((j, i), score)
            
    return scoringMatrix

# using the scoring matrix this method computes the optimal path matrix
def _computeOptimalPathMatrix(scoringMatrix, gapPentalty=-2):
    optimalScoreMatrix = np.zeros((len(scoringMatrix) + 1, len(scoringMatrix[0, ]) + 1))
    # matrix that represents the optimal route: 0 - diagonal, 1 - top, 2 - left
    optimalPathMatrix = np.zeros((len(scoringMatrix) + 1, len(scoringMatrix[0, ]) + 1))       
    
    for i, row in enumerate(optimalPathMatrix) :
        for j, elem in enumerate(optimalPathMatrix[i, ]) : 
            if i >= 1 :  # not in the first row anymore
                top = optimalScoreMatrix[i - 1, j] + gapPentalty
                if j >= 1 :  # not in the first row or column
                    topleft = optimalScoreMatrix[i - 1, j - 1] + scoringMatrix[i - 1, j - 1]
                    left = optimalScoreMatrix[i, j - 1] + gapPentalty
                    if topleft >= top and topleft >= left :
                        optimalScoreMatrix[i, j] = topleft
                        optimalPathMatrix[i, j] = 0
                    elif top >= topleft and top >= left :
                        optimalScoreMatrix[i, j] = top
                        optimalPathMatrix[i, j] = 1
                    else : 
                        optimalScoreMatrix[i, j] = left
                        optimalPathMatrix[i, j] = 2
                else :  # in the first column
                    optimalScoreMatrix[i, j] = top
                    optimalPathMatrix[i, j] = 1
            else : 
                if j >= 1 :  # in the first row, not first column
                    left = optimalScoreMatrix[i, j - 1] + gapPentalty
                    optimalScoreMatrix[i, j] = left
                    optimalPathMatrix[i, j] = 2
                else :  # in the top left
                    optimalScoreMatrix[i, j] = 0
                    optimalPathMatrix[i, j] = -1
                    
    print("optimal Scoring matrix" + str(optimalScoreMatrix))
    
    return optimalPathMatrix
       
# this method backTracks from the optimalPath matrix to get an alignment
# which is a list of the two sequences              
def _backTrackForAlignment(optimalPathMatrix, protein1, protein2):         
    alignment1 = ""
    alignment2 = ""
        
    currentScore = 0
    currentX = len(optimalPathMatrix[0, ]) - 1
    currentY = len(optimalPathMatrix) - 1
    while (currentX > 0 and currentY > 0) :
        currentScore = optimalPathMatrix.item((currentY, currentX))
        if (currentScore == 0) :  # diagonal, add both letters to alignment
            alignment1 = str(protein1[currentX - 1]) + alignment1
            alignment2 = str(protein2[currentY - 1]) + alignment2
            currentX -= 1
            currentY -= 1
        elif (currentScore == 1) :  # top, add only protein2 letter gap protein1
            alignment1 = "_" + alignment1
            alignment2 = str(protein2[currentY - 1]) + alignment2
            currentY -= 1       
        elif (currentScore == 2):  # left, add the top, gap the left
            alignment1 = str(protein1[currentX - 1]) + alignment1
            alignment2 = "_" + alignment2
            currentX -= 1
        else :  # at the end so stop
            break
    
    # add the traling letters and gaps
    if (currentX > 0) :
        while (currentX > 0) :
            alignment1 = str(protein1[currentX - 1]) + alignment1
            alignment2 = "_" + alignment2
            currentX -= 1
    elif (currentY > 0) :
        while (currentY > 0) :
            alignment2 = str(protein2[currentY - 1]) + alignment2
            alignment1 = "_" + alignment1
            currentY -= 1
    
    
    return [alignment1, alignment2]
            
def computeAlignmentScore(alignment, scoringMatrix, gapPentalty=-2,):    
    seq1 = alignment[0]
    seq2 = alignment[1]
    
    total_score = 0.0
    for i, letter1 in enumerate(seq1) :
        letter2 = seq2[i]
        if letter1 == '_' or letter2 == '_' :
            total_score += gapPentalty
        else : 
            total_score += scoringMatrix.get_value(letter1, letter2)
    
    return total_score
     
interactiveInputAlignment()

