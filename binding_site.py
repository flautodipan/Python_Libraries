"""
dato un frame, stampo la binding site
come la vuole gromacs per fare l'index
"""
#%%
import BioAlessandria as BA
import numpy as np
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

now_path    =   '../MD/WTC1/zerk/'
filename    =   'Protein_BS_WTC1.txt'
treshold    =   10 #angstrom