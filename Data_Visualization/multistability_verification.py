#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sit
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets
from IPython.display import display
from mpl_toolkits.mplot3d import Axes3D
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import pandas as pd
plt.rc("figure",figsize=(20,10))


# In[2]:


#Parametres of  the neutral
N,m,kmean,R0= 3,6.,0.5,5.
beta=m*R0
#elements clefs du neutral model
Ss,Ts,Is=1/R0,1-1/R0,(R0-1)/(R0+R0*kmean*(R0-1))
IIs,IIijs=Ts-Is,kmean*Is*Ts /Ss
N=3;
A0=np.random.random((N,N))
A0[0,0]=-0.;A0[0,1]=-0.5;A0[0,2]=+0.5
A0[1,0]=-0.5;A0[1,1]=+0.1;A0[1,2]=+0.5
A0[2,0]=+0.6;A0[2,1]=-0.5;A0[2,2]=+0.
def slowdyn(z,t,mu=Is/IIs,D=1,A=-A0):
    Nonsym=np.dot(A-A.T,z)
    sym=np.dot(A,z)
    Q=np.dot(z.T,np.dot(A,z))
    return D*z*((mu*Nonsym+sym-Q))
def slowdynr(z,t,mu=Is/IIs,D=1,A=-A0):
    Nonsym=np.dot(A-A.T,z)
    sym=np.dot(A,z)
    Q=np.dot(z.T,np.dot(A,z))
    return -D*z*((mu*Nonsym+sym-Q))


# In[3]:


##On commence donc par calculer ces solutions
def steady_state(A,mu):
    n=len(A)
    mat=mu*(A.T-A)+A.T
    C=np.concatenate( (mat[0:n-1,:]-mat[1:n,:],np.ones((1,n))), axis=0)
    det=np.linalg.det(C)
    if det==0: 
        print('Encore du travail !')
        return None 
    b=np.zeros((n,1))
    b[n-1]=1
    Xsol=np.linalg.solve(C,b)
    if len(Xsol[Xsol<0])>=1:
        print('Pas de solution positive')
        return None
    else : 
        print('Solution de coexistence')
        return Xsol
#Version alternative avec les fitness.
#Fonction de calcul des fitness
def fitness(A,mu):
    n=len(A)
    return mu*(A.T-A)+A.T-np.ones((n,1))@np.reshape(np.diagonal(A),(1,n))
##Calcul des solutions avec la version fitness
def steady_state_fit(FIT):
    n=len(FIT)        
    C=np.concatenate( (np.ones((1,n)),(FIT[0:n-1,:]-FIT[1:n,:])), axis=0)
    det=np.linalg.det(C)
    if det==0: 
        #print('Encore du travail !')
        return 0
    b=np.zeros((n,1))
    b[0]=1
    Xsol=np.linalg.solve(C,b)
    if len(Xsol[Xsol<0])>=1:
        #print('Pas de solution positive')
        return 0
    else : 
        #print('Solution de coexistence')
        return Xsol


# In[4]:


def Jac(FIT,z):
    n=len(FIT)
    z=np.reshape(z,(n,1))
    JAC_coex=np.diag(z.T[0])@(FIT-np.ones((n,n))@np.diag(np.array((FIT+FIT.T)@z).T[0]))
    bloc=np.array((FIT)@z -((z.T)@FIT@z)).T[0]
    JAC_trivial=np.diag(bloc)
    return JAC_coex+JAC_trivial


# In[5]:


def Jac_red(FIT,z): 
    n=len(FIT)
    s2,s3,s6=np.sqrt(2),np.sqrt(3),np.sqrt(6)
    P=np.concatenate( (np.ones((1,n)),
                       np.concatenate((np.ones((n-1,1)),
                                       -np.eye(n-1)),
                                      axis=1)),
                     axis=0)
    JAC=Jac(FIT,z)
    JACnew=np.linalg.inv(P)@JAC@P
    return JACnew[1:,1:]
def stab(FIT,z):
    M=Jac_red(FIT,z)
    VP=np.linalg.eig(M)[0]
    U=-1
    for v in VP:
        if v>0: U=1
        if (v==0&U==-1):U=0
    return U
def stab_val(FIT,z):
    M=Jac_red(FIT,z)
    VP=np.linalg.eig(M)[0]
    return np.max(VP.real)


# In[6]:


from itertools import combinations


# In[7]:


##Calculs des solutions stationnaires et de leurs stabilités
def classification_test(A=0,mu=0,FIT=[]):
    if len(FIT)==0:
        FIT=fitness(A,mu)
    N=len(FIT)
    ZZ,UU,ZZtype=[],[],[]


    PI=partiesliste(range(N))
    
    for I in PI:
        FITsub=np.zeros((len(I),len(I)))
        for i in range(len(I)):
            for j in range(len(I)):
                FITsub[i,j]=FIT[I[i],I[j]]
        zsub=steady_state_fit(FITsub)
        if not(isinstance(zsub, int)):
            z=np.zeros((N,1))
            z[I]=zsub
            UU.append(stab(FIT,np.array(z)))
            ZZ.append(z.reshape(N).tolist())
            ZZtype.append(len(I))
    return UU,ZZ,ZZtype
        
##Calculs des solutions stationnaires et de leurs stabilités
def classification_test_bis(A=0,mu=0,FIT=[]):
    if len(FIT)==0:
        FIT=fitness(A,mu)
    N=len(FIT)
    ZZ,UU,ZZtype,ZZeigen=[],[],[],[]


    PI=partiesliste(range(N))
    
    for I in PI:
        FITsub=np.zeros((len(I),len(I)))
        for i in range(len(I)):
            for j in range(len(I)):
                FITsub[i,j]=FIT[I[i],I[j]]
        zsub=steady_state_fit(FITsub)
        if not(isinstance(zsub, int)):
            z=np.zeros((N,1))
            z[I]=zsub
            UU.append(stab(FIT,np.array(z)))
            ZZ.append(z.reshape(N).tolist())
            ZZtype.append(len(I))
            ZZeigen.append(stab_val(FIT,np.array(z)))
    return UU,ZZ,ZZtype,ZZeigen


# In[8]:


def classification_test_print(FIT,s):
    N=len(FIT)
    ZZ=[]
    for I in combinations(range(N),s):
        FITsub=np.zeros((len(I),len(I)))
        for i in range(len(I)):
            for j in range(len(I)):
                FITsub[i,j]=FIT[I[i],I[j]]
        zsub=steady_state_fit(FITsub)
        if not(isinstance(zsub, int)):
            z=np.zeros((N,1))
            z[list(I)]=zsub
            u=stab(FIT,np.array(z))
            if u<0 :
                ZZ.append(z.reshape(N).tolist())
    return ZZ


# In[11]:


def resul_gen(FIT, s0=1,s1=2, fprint=False):
    txt=''
    for s in range(s0,s1+1):
        ZZ=classification_test_print(FIT,s)
        p=len(ZZ)
        if p>1 : 
            chain=f'Il y a {p} états stables de {s} souche'+'s'*(s>1)+' : \n'
            for z in ZZ :
                chain+=str(z)+'\n'
        elif p>0 : 
            chain=f'Il y a {p} unique état stable de {s} souche'+'s'*(s>1)+' : \n'
            chain+=str(ZZ[0])+'\n'
        else : chain=f'Il n\'y a pas d\'état stable de {s} souche'+'s'*(s>1)
        if fprint : print(chain)
        txt+=chain+'\n'
        
    return txt

