'''
Created on Apr 24, 2017

@author: Goldacbj
'''
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import scipy.stats as stats

def lesson424():
    # mutlitple alignement and droping gaps for one sequence alignment
    # then creating a position specific mutation matrix 
    align1 = AlignIO.read('file', "fasta")
    print(str(align1[0:3]))
    
    id = (protein.id for protein in align1)
    
    align_array = np.array([list(protein.seq) for protein in align1], dtype='str', order='F')
    
    align2 = pd.DataFrame(dtat=align_array, index=id)
    align2.iloc[0:3, 200:300]
    align2.iloc[:, 205] # coulmn 205
    
    # remove all the gaps from human 
    (m, n) = align2.shape
    # n should get down to 110
    
    drop_list = []
    for k in range(0, n) :
        if align2.iloc[0, k] == '-' : 
            drop_list.append(k)
            
    align3 = align2.drop(drop_list, axis=1) # human align without any gaps
    
    align3.iloc[:, 0] # the first letter of each 
    
    ### create a specific position mutation matrix 
    align3[167].value_counts() # counts for each letter into a map of letter and count 
    ## get a mutation matrix from these counts 
    # should drop first row, since we are comparing to this 
    human_seq = align3.iloc[1, :]
    align = align3.drop(align3.index[0], axis=0) #drop row
    
    counts = align.apply(pd.value_counts).fillna(0).T # functional prog # secquence down left side with each letter count in a column 
    
    # mutMatrix.ilc[0, :] # mutation rate of each amino acid 
    
    #mutability = np.sezors(len(human_seq))


    # show in a bar graph ---------------- 
    # mutatMatix.iloc[100::-1, "].plot.harh(stacke=True); plt.show()
    
    # compute mutabilit of each ammino acid ----------------------------- 
    # mutmat = counts.div(counts.sum(axis=1), axis=0)
    # mutabilitiy[k] = 1 - mutatMatrix.iloc[k, :][human_seeq.values[2]]
    # mutablilit = pd.datafram(mutabiliti, index= human_seq.values)
    # mutability .plot.bar(); plt.show();
    
    # compute position specific scoring matrix
    
    # compute distribution matrix by alignment column 
    #    aaprob = counts.sum()/counts.sum().sum()
    #    logs odd ratio of aminoAcid at position k vs just a random acid at position k  (multiply by 10 to get scoring matrix)
    #    AA[insulin] = mutmat.iloc[k, :]
    #    AA['random'] = aaprobs
    #     aa.plot.bar(); plt.show()
    
    # aA.iloc[: 0]
    # PSM = np.ound(10 * np.log10(mutmat.div(AAprob.iloc[:, 0], axis=1)))
    
    
    
    
    