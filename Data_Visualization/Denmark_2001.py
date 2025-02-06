#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# In[2]:
data = pd.read_csv('Denmark_Identities_2001.csv')


# In[2]:


z_Denmark = np.array([92, 86, 53, 45, 27,  22, 21, 19, 18, 14, 12, 10, 10,
                       8, 8, 7, 5, 4, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1])
z_Denmark = z_Denmark/sum(z_Denmark)


# In[3]:


Denmark_2001 = {'strain': data['Denmark'],
                'cases':  z_Denmark}
Denmark_2001 = pd.DataFrame(Denmark_2001).sort_values(['cases'], ascending = False)
Denmark_2001.strain.astype('str')


# In[4]:


def D2001():
    return Denmark_2001

