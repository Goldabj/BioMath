'''
Created on May 14, 2017

@author: Goldacbj
'''
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import scipy.stats as stats
from multiprocessing import pool as pl

def main():
    # ----------------------- PART 1 ---------------------------
    ### carrier Proteomes ----------------------
    carrier_proteomes = dict()
    
    # Anopheles-arabiensis-Dongola_PEPTIDES_AaraD1.3.fa
    carrier1_proteome =  list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-arabiensis-Dongola_PEPTIDES_AaraD1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-arabiensis-Dongola_PEPTIDES_AaraD1.3.fa'] = carrier1_proteome

    # Anopheles-darlingi-Coari_PEPTIDES_AdarC3.3.fa
    carrier2_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-darlingi-Coari_PEPTIDES_AdarC3.3.fa", "fasta"))
    carrier_proteomes['Anopheles-darlingi-Coari_PEPTIDES_AdarC3.3.fa'] = carrier2_proteome

    # Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa
    carrier3_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa", "fasta"))
    carrier_proteomes['Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa'] = carrier3_proteome

    # Anopheles-atroparvus-EBRO_PEPTIDES_AatrE1.3.fa
    carrier4_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-atroparvus-EBRO_PEPTIDES_AatrE1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-atroparvus-EBRO_PEPTIDES_AatrE1.3.fa'] = carrier4_proteome

    # Anopheles-coluzzii-Mali-NIH_PEPTIDES_AcolM1.2.fa
    carrier5_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-coluzzii-Mali-NIH_PEPTIDES_AcolM1.2.fa", "fasta"))
    carrier_proteomes['Anopheles-coluzzii-Mali-NIH_PEPTIDES_AcolM1.2.fa'] = carrier5_proteome

    # Anopheles-culicifacies-A-37_PEPTIDES_AculA1.3.fa
    carrier6_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-culicifacies-A-37_PEPTIDES_AculA1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-culicifacies-A-37_PEPTIDES_AculA1.3.fa'] = carrier6_proteome

    # Anopheles-dirus-WRAIR2_PEPTIDES_AdirW1.3.fa
    carrier7_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-dirus-WRAIR2_PEPTIDES_AdirW1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-dirus-WRAIR2_PEPTIDES_AdirW1.3.fa'] = carrier7_proteome

    # Anopheles-farauti-FAR1_PEPTIDES_AfarF2.1.fa
    carrier8_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-farauti-FAR1_PEPTIDES_AfarF2.1.fa", "fasta"))
    carrier_proteomes['Anopheles-farauti-FAR1_PEPTIDES_AfarF2.1.fa'] = carrier8_proteome

    # Anopheles-funestus-FUMOZ_PEPTIDES_AfunF1.3.fa
    carrier9_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-funestus-FUMOZ_PEPTIDES_AfunF1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-funestus-FUMOZ_PEPTIDES_AfunF1.3.fa'] = carrier9_proteome

    # Anopheles-maculatus-maculatus3_PEPTIDES_AmacM1.3.fa
    carrier10_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-maculatus-maculatus3_PEPTIDES_AmacM1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-maculatus-maculatus3_PEPTIDES_AmacM1.3.fa'] = carrier10_proteome

    # Anopheles-sinensis-SINENSIS_PEPTIDES_AsinS2.1.fa
    carrier11_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-sinensis-SINENSIS_PEPTIDES_AsinS2.1.fa", "fasta"))
    carrier_proteomes['Anopheles-sinensis-SINENSIS_PEPTIDES_AsinS2.1.fa'] = carrier11_proteome

    # Anopheles-stephensi_SDA-500_PEPTIDES.AsteS1.3.fa
    carrier12_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-stephensi_SDA-500_PEPTIDES.AsteS1.3.fa", "fasta"))
    carrier_proteomes['Anopheles-stephensi_SDA-500_PEPTIDES.AsteS1.3.fa'] = carrier12_proteome

    ### Mild carrier Proteomes -----------------------------------------
    mild_carrier_proteomes = dict()
    
    # Anopheles-albimanus_STECLA_PEPTIDES_AalbS1.3.fa
    mild_carrier1_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-albimanus_STECLA_PEPTIDES_AalbS1.3.fa", "fasta"))
    mild_carrier_proteomes['Anopheles-albimanus_STECLA_PEPTIDES_AalbS1.3.fa'] = mild_carrier1_proteome

    # Anopheles-epiroticus-Epiroticus2_PEPTIDES_AepiE1.3.fa
    mild_carrier2_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-epiroticus-Epiroticus2_PEPTIDES_AepiE1.3.fa", "fasta"))
    mild_carrier_proteomes['Anopheles-epiroticus-Epiroticus2_PEPTIDES_AepiE1.3.fa'] = mild_carrier2_proteome

    # Anopheles-melas-CM1001059_PEPTIDES_AmelC2.1.fa
    mild_carrier3_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-melas-CM1001059_PEPTIDES_AmelC2.1.fa", "fasta"))
    mild_carrier_proteomes['Anopheles-melas-CM1001059_PEPTIDES_AmelC2.1.fa'] = mild_carrier3_proteome

    # Anopheles-merus-MAF_PEPTIDES_AmerM2.1.fa
    mild_carrier4_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-merus-MAF_PEPTIDES_AmerM2.1.fa", "fasta"))
    mild_carrier_proteomes['Anopheles-merus-MAF_PEPTIDES_AmerM2.1.fa'] = mild_carrier4_proteome

    # Anopheles-minimus-Minimus1_PEPTIDES_AminM1.3.fa
    mild_carrier5_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-minimus-Minimus1_PEPTIDES_AminM1.3.fa", "fasta"))
    mild_carrier_proteomes['Anopheles-minimus-Minimus1_PEPTIDES_AminM1.3.fa'] = mild_carrier5_proteome
    
    ### Non Carrier Proteomes ------------------
    non_carrier_proteomes = dict()
    
    # Anopheles-quadriannulatus-SANGQUA_PEPTIDES.AquaS1.3.fa
    non_carrier1_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-quadriannulatus-SANGQUA_PEPTIDES.AquaS1.3.fa", "fasta"))
    non_carrier_proteomes['Anopheles-quadriannulatus-SANGQUA_PEPTIDES.AquaS1.3.fa'] = non_carrier1_proteome

    # Anopheles-christyi-ACHKN1017_PEPTIDES_AchrA1.3.fa
    non_carrier2_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/Proteomes/Anopheles-christyi-ACHKN1017_PEPTIDES_AchrA1.3.fa", "fasta"))
    non_carrier_proteomes['Anopheles-christyi-ACHKN1017_PEPTIDES_AchrA1.3.fa'] = non_carrier2_proteome
    
    
    # Proteins to compare -----------------------     
    protein_sequences = list()

    # Protein 1: 
    # lysysome protein 
    #    # MKVFIAIVLTIVASCALAEAKKFSKCDLAKTLANNGIARASL
    #    # PDWICLVQNESAFSTSATNKNKNGSTDYGIFQINNKYWCDSSYGSNDCK
    #    # IACKKLLDDDITDDIKCAKLIFKRHGYNAWYGWKNHCNGKALPNVDSCF
    protein1 = SeqRecord(Seq("MKVFIAIVLTIVASCALAEAKKFSKCDLAKTLANNGIARASLPDWICLVQNESAFSTSATNKNKNGST" 
        + "DYGIFQINNKYWCDSSYGSNDCKIACKKLLDDDITDDIKCAKLIFKRHGYNAWYGWKNHCNGKALPNVDSCF",  IUPAC.IUPACProtein), id="ApLysc1")
    protein_sequences.append(protein1)
    # Protein 2: 
    # F2Y0G7
    protein2 = SeqIO.read("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/F2Y0G7.fasta.txt", "fasta")
    protein_sequences.append(protein2)
    # Protein 3: 
    # F5HM55
    protein3 = SeqIO.read("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/F5HM55.fasta.txt", "fasta")
    protein_sequences.append(protein3)
    # Protein 4: 
    # L7RM09
    protein4 = SeqIO.read("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/Project/L7RM09.fasta.txt", "fasta")
    protein_sequences.append(protein4)
    
    outputFile = open("output.txt", "r+")
    
    s = computeOptimalAlignmentScores(protein_sequences, outputFile)
    print(str(s))
    
    results = list()
    
    results = findMaxAlignmentForProteins(protein_sequences, non_carrier_proteomes, "Non-Carrier", 8)
    writeAlignmentResults(results, outputFile)
    print(str(results))
    
    results = findMaxAlignmentForProteins(protein_sequences, mild_carrier_proteomes, "Mild-Carrier", 8)
    writeAlignmentResults(results, outputFile)
    print(str(results))
    
    results = findMaxAlignmentForProteins(protein_sequences, carrier_proteomes, "Carrier", 8)
    writeAlignmentResults(results, outputFile)
    print(str(results))
    
    outputFile.close()



def computeOptimalScore(protein_seq, scoring_matrix):
    '''Computes the max alignemtns score of a protein sequence'''
    score = 0
    for letter in protein_seq : 
        score = score + scoring_matrix[letter, letter]
    
    return score
        
        
def computeOptimalAlignmentScores(protein_records, output_file):
    '''writes the optimal alignment Scores to the ouptputfile'''
    final_string = ""
    scoringMatrix = matlist.blosum62
    output_file.write("optimal scores: --------------------------------- \n")
    for record in protein_records : 
        score = computeOptimalScore(record.seq, scoringMatrix)
        output_file.write("\t" + record.id + " : " + str(score) + "\n")
        final_string = final_string + "\t" + record.id + " : " + str(score) + "\n"
    
    output_file.write("\n + \n")
    return final_string


def computeMaxAlignment(protein1, proteome, proteome_id, proteome_type):
    best_alignment = ("1", "2", 3, 0, 0)
    highest_score = -10000
    for protein in proteome: 
        try: 
            alignment_score = pairwise2.align.globalds(str(protein1.seq), str(protein.seq), matlist.blosum62, -11, -1, score_only=True)
            if alignment_score > highest_score : 
                highest_score = alignment_score
                best_alignment = pairwise2.align.globalds(str(protein1.seq), str(protein.seq), matlist.blosum62, -11, -1, one_alignment_only=True)
        except Exception as inst:
            print(inst)
            print("\n")
    
    return (proteome_type, proteome_id, best_alignment)

def writeAlignmentResults(results, output_file) :
    # TODO: format the output in the file better - by type of proteome 
    output_file.write("\n")
    output_file.write("optimal alignments ---------------------------------- \n")
    for result in results: 
        output_file.write("\t \t" + str(result))
        output_file.write("\n")

    return 

def findMaxAlignmentForProteins(protein_records, proteomes, proteomes_type, numofThreads):    
    pool = pl.Pool(processes=numofThreads)
    tasks = list()
    
    # tasks should be each proteome for each protein
    for protein in protein_records : 
        for key in proteomes : 
            tasks.append((protein, proteomes[key], key, proteomes_type))
    
    results = list()
    for task in tasks :
        results.append(pool.apply_async(computeMaxAlignment, task))
    
    pool.close()
    pool.join()
    
    return results

