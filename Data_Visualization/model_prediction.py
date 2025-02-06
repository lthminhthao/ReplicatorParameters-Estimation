#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy import integrate


# In[2]:


#FUNCTION SOLVES REPLICATOR SYSTEM
def RE(n, lamd):   
    def du_dt(u, t = 0):
        return u*(np.dot(lamd,u) - np.dot(np.dot(lamd,u),u))
               
    t = np.linspace(0, 50, 500)
    #u0 = np.random.random(n)
    #u0 = u0/sum(u0)
    u0 = (1/n)*np.ones(n)
    #solve and compute the time
    u, infodict = integrate.odeint(du_dt, u0, t, full_output=True)
    infodict['message']
    #solution
    sol = u.T[:,-1]
    return(sol)


# In[3]:


#Solve the replicator system with u0 as a parameter
def Replicator(n, lamd, time, u0):
    def du_dt(u, t = 0):
        return u*(np.dot(lamd,u) - np.dot(np.dot(lamd,u),u))
               
    t = np.linspace(0, time, 500)
    u, infodict = integrate.odeint(du_dt, u0, t, full_output=True)
    infodict['message']
    #solution
    sol = u.T[:,-1]
    return(sol)


# In[4]:


def matrix(n, create, mean, var):
    list_strain = [] #list of fequencies for 70 strains
    for i in range (0, n):
        list_strain.append([])
        
    min_matrix = []
    max_matrix = []
    mean_matrix = []
    
    
    for test in range (0, 50):
    ### Create a symmetric normally-distributed random fitness matrix for each country
        lambda_matrix = create(n, mean, var)
    #Solve Replicator systems for each country    
        solution = RE(n, lambda_matrix)
        solution = sorted(solution, reverse=True) #solution for a country in one test
        
        for j in range (0, n): #list of fequencies for 70 strains saved replicator solution for all tests for one country
            list_strain[j].append(solution[j]) #list has 70 strains, each elements contains all serotypes frequencies of that strain through all tests
    
    for i in range (0, n):
        min_matrix.append(np.min(list_strain[i]))
        mean_matrix.append(np.mean(list_strain[i]))
        max_matrix.append(np.max(list_strain[i]))
        
    matrix = np.array([min_matrix, mean_matrix, max_matrix])
    
    return matrix

