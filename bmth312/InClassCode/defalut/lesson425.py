'''
Created on Apr 25, 2017

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


def lesson425():
    # sequence logos 
    