'''
Created on Mar 30, 2017

@author: Goldacbj
'''
import threading
import time
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist 

def main() : 
    
    human_insulin = SeqIO.read("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/insulin_human.fasta.txt", "fasta")
    human_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/proteome_human.fasta.txt", "fasta"))
    fruityfly_proteome = list(SeqIO.parse("C:/Users/Goldacbj/Google Drive/Documents/BIO/BMTH312/proteome_fruitfly.fasta.txt", "fasta")) 
    
    startTime = time.time()
    ortholog = findOrtholog(protein=human_insulin, proteomeLookIn=fruityfly_proteome, proteomeFrom=human_proteome)
    endtime = time.time()
    totaltime = endtime - startTime
    print("Search time = " + str(totaltime))
    
    if ortholog == None: 
        print("Ortholog found: " + str(ortholog))
    else: 
        print("Ortholog not found")
        
def _searchForBestMatch(sequence, proteome, startIndex, finishIndex, possibleCanadates):
    ''' returns the max sequence match, the score, and index in proteome, looking in 
    the proteome form the startIndex to the finish Index'''
    print("Thread started searching through " + str(startIndex) + " -- " + str(finishIndex))
    SubMat = matlist.blosum62 
    maxIndex = 0
    maxScore = 0
    score = 0
    for i in range(startIndex, finishIndex) : 
        seq2 = proteome[i].seq
        if (i % 500 == 0) :
            print("processing#" + str(i) + "...")
        try: 
            score = pairwise2.align.localds(sequence, seq2, SubMat, -11, -1, score_only=True)
        except SystemError : 
            print("error found: " + str(i))
            continue
        if score > maxScore : 
            maxScore = score
            maxIndex = i
    
    possibleCanadates.append((proteome[maxIndex].seq, maxScore, maxIndex))
    print("Thread Finsihed searching through " + str(startIndex) + " -- " + str(finishIndex))
    

def searchForBestMatchThreaded(sequence, proteome, numOfThreads):
    ''' search for best match from the proteome using 4 threads. returns list with (sequence, score, indexOfProtein)'''
    possibleCanadates = list()
    threads = list()
    split = round(len(proteome) / numOfThreads)
    for i in range(0, numOfThreads) :
        thread = threading.Thread(target=_searchForBestMatch, args=(sequence, proteome, split * i, split * (i + 1), possibleCanadates))
        threads.append(thread)
        thread.start()
    
    for thread in threads : 
        thread.join()
    
    bestCanadate = ()
    maxScore = 0
    for possible in possibleCanadates : 
        possibleScore = possible[1]
        if (possibleScore > maxScore) :
            maxScore = possibleScore
            bestCanadate = possible
            
    return bestCanadate


def findOrtholog(protein, proteomeLookIn, proteomeFrom):
    '''finds the ortholog from the first protein in the proteome
    which we are searching in using the optimal alignment and 
    blosum 62 scoring matrix'''
    
     # # search for canadaites in proteomeLookIn 
    canadate = searchForBestMatchThreaded(protein.seq, proteomeLookIn, 6)    
    print("best canadate: " + str(canadate))
    
    # # search for match to canadate in protemoFrom
    bestMatch = searchForBestMatchThreaded(canadate[0], proteomeFrom, 6)
    print("best match for canadate: " + str(bestMatch))
    
    # # check to see if match is the original protein
    if (str(bestMatch) == (str(protein.seq))) :
        return proteomeLookIn[canadate[2]]
    else : 
        return
    

main()
