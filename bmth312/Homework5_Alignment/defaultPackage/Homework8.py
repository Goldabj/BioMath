'''
Created on Apr 3, 2017

@author: Goldacbj
'''
from Bio import SeqIO
from Bio import pairwise2
import numpy as np
import pandas as pd
from Bio.SubsMat import MatrixInfo as matlist 

def main():
    human_insulin = SeqIO.read("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/insulin_human.fasta.txt", "fasta")
    other_insulins = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/insulin_101_uniprot.fasta.txt", "fasta"))
    blosum80 = matlist.blosum80;
    
    # A computePseduoCountMatrix
    countMatrix = computePseudoCountMatrix(protein=human_insulin, otherSequences=other_insulins, subsMatrix=blosum80, pseudoCount=1)
    print("A: Pseudo Count Matrix: ")
    print(str(countMatrix))
    
    # B: computePusedoCountMatrix(pseudo=0), & computMutationRate
    print("Part B:  ")
    print(" ")
    countMatrix = computePseudoCountMatrix(protein=human_insulin, otherSequences=other_insulins, subsMatrix=blosum80, pseudoCount=0)
    mutationList = computeMutationRate(countMatrix)
    print("Mutation rates: ")
    print(str(mutationList))
    
    # C : for each identity, getSeqIdentiyAlginements and computeMutationRate
    print("Part C: ")
    print("")
    zeroPercentRecords = other_insulins
    twentyPercentRecords = getSeqIdentityAlignments(protein=human_insulin, otherSequences=other_insulins, subsMatrix=blosum80, sequenceIdentity=0.20)
    fiftyPercentRecords = getSeqIdentityAlignments(protein=human_insulin, otherSequences=other_insulins, subsMatrix=blosum80, sequenceIdentity=0.50)
    nintyPercentRecords = getSeqIdentityAlignments(protein=human_insulin, otherSequences=other_insulins, subsMatrix=blosum80, sequenceIdentity=0.90)
    
    totalcountMatrix = computePseudoCountMatrix(protein=human_insulin, otherSequences=zeroPercentRecords, subsMatrix=blosum80, pseudoCount=0)
    twentyCountMatrix = computePseudoCountMatrix(protein=human_insulin, otherSequences=twentyPercentRecords, subsMatrix=blosum80, pseudoCount=0)
    fiftyCountMatrix = computePseudoCountMatrix(protein=human_insulin, otherSequences=fiftyPercentRecords, subsMatrix=blosum80, pseudoCount=0)
    nintyCountMatrix = computePseudoCountMatrix(protein=human_insulin, otherSequences=nintyPercentRecords, subsMatrix=blosum80, pseudoCount=0)
    
    totalMutationRate = computeMutationRate(totalcountMatrix).sum()
    twentyMutationRate = computeMutationRate(twentyCountMatrix).sum()
    fiftyMutationRate = computeMutationRate(fiftyCountMatrix).sum()
    nintyMutationRate = computeMutationRate(nintyCountMatrix).sum()
    
    print("0% sequence Id mutationRate:  " + str(totalMutationRate))
    print("20% sequence Id mutationRate:  " + str(twentyMutationRate))
    print("50% sequence Id mutationRate:  " + str(fiftyMutationRate))
    print("90% sequence Id mutationRate:  " + str(nintyMutationRate))

def computePseudoCountMatrix(protein, otherSequences, subsMatrix, pseudoCount=1):
    ''' aligns all the other sequences with the sequence, then constructs a count matrix
    with the alignment matches. Then adds the psedeoCount to the count matrix and returns 
    that matrix as a pandas matrix'''
    aa_letters = list("ARNDCQEGHILKMFPSTWYV")
    countMatrix = pd.DataFrame(np.zeros((20,20)),index=aa_letters,columns=aa_letters)
    
    seq1 = protein.seq
    
    for k, record in enumerate(otherSequences) : 
        seq2 = record.seq
        try: 
            alignment = pairwise2.align.globalds(seq1, seq2, subsMatrix, -11, -1, one_alignment_only=True)        
            alignmentSeq1 = alignment[0][0]
            alignmentSeq2 = alignment[0][1]
    
            for i, letter in enumerate(alignmentSeq1): 
                letter2 = alignmentSeq2[i]
                count = countMatrix.get_value(letter, letter2)
                if (count != None) : 
                    countMatrix.set_value(letter, letter2, count+1)
        except: 
            'print("Error in aligning seq: " + str(k))'
        
    countMatrix = countMatrix + pseudoCount
    return countMatrix 
            
        
    
    
def computeMutationRate(countMatrix):
    ''' takes in a count matrix and then computes the mutablility for each 20 amino acid.
    it returns the mutation rates in order of lowest highest'''
    symmetricCountMatrix = 0.5 * (countMatrix + countMatrix.T) # symmetric counts
    probabilityMatrix = symmetricCountMatrix/symmetricCountMatrix.sum().sum() # P(A n B) 
    individualProbs = probabilityMatrix.sum() # individual probs -- P(A), P(B) etc
    mutationMatrix = probabilityMatrix.div(individualProbs, axis=1)
    
    mutationRates = pd.Series(1-np.diag(mutationMatrix), index=mutationMatrix.columns) 
    mutationRates = mutationRates * individualProbs
    
    return mutationRates.sort_values()
    
    
def getSeqIdentityAlignments(protein, otherSequences, subsMatrix, sequenceIdentity):
    ''' returns all the alignemtns which equals or are above the sequence Identity rate, returns 
    a list of the good records'''
    goodSequences = list()
    seq1 = protein.seq
    for k, record in enumerate(otherSequences) : 
        seq2 = record.seq

        alignment = pairwise2.align.globalds(seq1, seq2, subsMatrix, -11, -1, one_alignment_only=True)  
            
        seqIdentity = computeSeqIdentity(alignment[0])
        if (seqIdentity >= sequenceIdentity) : 
            goodSequences.append(record)
    
    return goodSequences

def computeSeqIdentity(alignment):
    seq1 = alignment[0]
    seq2 = alignment[1]
    
    minLength = min(len(seq1), len(seq2))
    matchCount = 0
    length = 0
    
    for i in range(0, minLength): 
        letter1 = seq1[i]
        letter2 = seq2[i]
        if (letter1 != "_") : 
            length = length + 1
            if (letter1 == letter2) : 
                matchCount = matchCount + 1
    
    return matchCount/length
    


main()