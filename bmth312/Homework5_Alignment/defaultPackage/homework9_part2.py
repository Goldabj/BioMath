'''
Created on Apr 7, 2017

@author: Goldacbj
'''

import defaultPackage.Homework9_statistics as part1
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import scipy.stats as stats
from multiprocessing import pool

def main():
    num2()
    

def num2():
    # searching for homologs and plotting their distribution 
    human_insulin = part1.getSequenceFromFile("sp|P01308|INS_HUMAN")
    fruitFlyProteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/proteome_fruitfly.fasta.txt", "fasta"))
    print("Searching for homologs")
    alignmentScores = computeAllAlignmentScores(human_insulin, fruitFlyProteome) # map of IDs and their alignment scores
    print("")
    plt.hist(list(alignmentScores.values()))
    plt.draw()
    
    # computing mean and std
    mean = np.mean(list(alignmentScores.values()))
    std = np.std(list(alignmentScores.values()))
    print("mean of scores: " + str(mean))
    print("std of scores: " + str(std))
    
    # z scores for each alignment
    z_scores = dict()
    for key in alignmentScores.keys() : 
        val = alignmentScores.get(key)
        z_scores[key] = (val - mean)/std
        
    print("z scores for alignments: " + str(list(z_scores.values()) ))
    print("num of z scores: " + str(len(list(z_scores.values()))))
    
    # p values 
    p_values = dict()
    for key in z_scores.keys() : 
        score = z_scores.get(key)
        p_values[key] = 1 - stats.norm.cdf(score)
    
    print("p values: " + str(p_values))
    
    
    # compute E -values and print possible homologs
    # computeing from len of p_values since some alignments didn't work .
    E_values = dict()
    for key in p_values.keys(): 
        p_val = p_values.get(key)
        E_values[key] = (len(p_values) * p_val)
        
    homologs = list()
    max = -1000;
    for key in E_values.keys() : 
        eVal = E_values.get(key)
        if (eVal > max) :
            homologs = list()
            max = eVal
            homologs.append(key)
        elif eVal == max : 
            homologs.append(key)
    
    
    # print homologs
    print("possible homologs: " + str(homologs))
    
    plt.show()
    

    

def computeAllAlignmentScores(seq, proteome):
    ''' returns all the sequneceIDs and their alignment scores in a dictionary'''
    scores = dict()
    tasks = list()
    scaler = round(len(proteome)/4)
    for i in range(0, 3): 
        tasks.append((seq, proteome, i*scaler, (i+1)*scaler))
    tasks.append((seq, proteome, 3*scaler, len(proteome)))
    
    poool = pool.Pool(processes=4)
    for i in range(0, 4) : 
        scores.update(poool.apply_async(computeAlignmentScore, tasks[i]).get())
        
    poool.close()
    poool.join()
    
    return scores

def computeAlignmentScore(seq, proteome, startIndex, endIndex):
    ''' computes alignment scores and returns a map of the sequence ID and the score'''
    scores = dict()
    subsMatrix = matlist.blosum62
    for i in range(startIndex, endIndex) : 
        seq2 = proteome[i].seq
        id = proteome[i].id
        try :
            score =  pairwise2.align.localds(seq, seq2, subsMatrix, -11, -1, score_only=True)
        except: 
            'nothing'
        scores[id] = score
    
    return scores

if __name__ == '__main__':
    main()