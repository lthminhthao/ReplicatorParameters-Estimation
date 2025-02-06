#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


def jacobian_func(fit, z):
    #input: fitness matrix, time, z(t)

    n = len(z)
    jacobian = np.zeros((n, n))
        
    for i in range(0, n):
        for j in range (0, n):
            
            #non-diagonal elements
            if j != i:
                jacobian[i, j] = fit[i,j]*z[i] - z[i]*np.matmul((fit.T + fit)[j, :], z)
            
            #diagonal elements
            if j == i:
                jacobian[i, i] = np.matmul(fit[i, :], z) - np.matmul(np.matmul(fit, z), z) - z[i]*np.matmul((fit.T + fit)[i, :], z)
    
    return jacobian

