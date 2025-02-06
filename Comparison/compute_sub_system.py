#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from model_prediction import Replicator


# In[2]:


#import list of all strains
all_strain = ['16', '23B', '33B', '7C', '23F', 'NT3b', '9V', '5', '35F', '24F', '12F', '9', '9L', 'NT2', '20', '45', '35A', 
              '19C', '10', '24A', '19F', '33C', 'NT4a', '40', '15A', '6D', '7B', '6A', '4', '11A', '33A', '46', '4B', '23A', 
              '19', '18A', '6A/6B', '11B', '38', '33F', '29', '18F', '13', '15C', 'NT', '6C', '18B', '11D', '17A', '21', '31', 
              '22F', '36', '8', '9A', '19B', '28F', '24B', '48', '6B', '23C', '10A', '12F/12A/4', '25F', '35B', '19A', '18C', 
              '39', '17F', 'NT4b', '28A', '1', '22A', '32F', '6', '35C', '15B', '37', '17', '18', 'NC', '15', '23', '11', 
              '10F', '14', '3', '9N', '16F', '34', '22', '7F']

#import big_alpha
big_alpha = np.loadtxt(open("big_alpha.txt","rb"), delimiter = ",", skiprows = 0)


# In[3]:


# With input a list of strain names (subset of big N), it extracts the relevant alpha submatrix from the big alpha matrix
# function extract alpha matrix for a given set of strains

def alpha_extract(strain_list):
    s = len(strain_list)
    alpha = np.zeros((s,s))

    for i in range (0, s):
        invader = strain_list[i]
        row = all_strain.index(invader)
        
        for j in range (0, s):
            resident = strain_list[j]
            col = all_strain.index(resident)
            
            alpha[i,j] = big_alpha[row, col]
    return alpha

