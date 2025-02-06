#!/usr/bin/env python
# coding: utf-8

# In[1]:


##Packages##
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib as mpl

##figures inline##
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rc("figure",figsize=(16,10))


# In[2]:


def network(LL):
    """This use the package networkx
    return the graph G from a given fitnesses matrix LL"""
    if LL is list :
        LL=np.array(LL)
    LL=LL-np.diag(np.diag(LL)) # be sure that the self fitness is zero
    G=nx.DiGraph()
    n=len(LL)
    nodes=['{}'.format(i+1) for i in range(n)]
    G.add_nodes_from(nodes)
    pair=[(i,j) for j in range(n) for i in range(j) 
           if LL[i,j] !=0 or LL[j,i] !=0  ]
    edges=[('{}'.format(x[0]+1),'{}'.format(x[1]+1)) for x in pair]
    G.add_edges_from(edges)
    for (i,j),(a,b) in zip(pair,G.edges):
        classe='degenerate'
        if LL[i,j]<0 and LL[j,i]<0: G[a][b]['type']='bist'
        elif LL[i,j]>0 and LL[j,i]>0: G[a][b]['type']='coex'
        elif LL[i,j]<0 and not LL[j,i]<0: G[a][b]['type']='excl_'+b
        elif not LL[i,j]<0 and LL[j,i]<0: G[a][b]['type']='excl_'+a
        #print(classe)
    return G


# In[3]:


def net_draw(G,**kwargs):
    """plot the graph G
    kwargs as the following entry :
    - nodes : z, the size of the nodes. Default is 100
    -ax= plot on the curent ax. default is None"""
    
    ax=kwargs.pop('ax', None) #name of axis
    if ax is None : fig,ax=plt.subplots()
    
    zz=np.array(kwargs.pop('nodes', [1/G.number_of_nodes()]*G.number_of_nodes()))#size of nodes
    
    edges_size=kwargs.pop('edges_size', None)#size of edges
    
    fs=kwargs.pop('font_size', 24)#size of font
    
    #topology of the graph
    pos=nx.layout.circular_layout(G)
    nx.draw_networkx(G,pos=pos,
            node_size=2000*zz+0.01,node_color='gainsboro', alpha = 1,
           edgelist=[], ax=ax,font_size=fs)
    ##definitions of the edges by types
    edgecoex,edgebist,edgeexclui,edgeexcluj=[],[],[],[]
    edgecoex_size,edgebist_size,edgeexclui_size,edgeexcluj_size=[],[],[],[]
    for i,j in G.edges :
        
        if G[i][j]['type']=='coex':
            edgecoex.append((i,j))
            edgecoex_size.append(20*zz[eval(i)-1]*zz[eval(j)-1])
        elif G[i][j]['type']=='bist':
            edgebist.append((i,j))
            edgebist_size.append(20*zz[eval(i)-1]*zz[eval(j)-1])
        elif G[i][j]['type']=='excl_'+i:
            edgeexclui.append((i,j))
            edgeexclui_size.append(20*zz[eval(i)-1]*zz[eval(j)-1])
        elif G[i][j]['type']=='excl_'+j:
            edgeexcluj.append((j,i))
            edgeexcluj_size.append(20*zz[eval(i)-1]*zz[eval(j)-1])
    ##size of the edges
    if edges_size=='ztype' : print('adaptative edges width')
    elif edges_size==None : 
        edgecoex_size=5
        edgebist_size=5
        edgeexclui_size=[5]*len(edgeexclui)
        edgeexcluj_size=[5]*len(edgeexcluj)
    else: 
        edgecoex_size=edges_size
        edgebist_size=edges_size
        edgeexclui_size=[edges_size]*len(edgeexclui)
        edgeexcluj_size=[edges_size]*len(edgeexcluj)
    ###final ploting###
    ##coexistence edges##
    nx.draw_networkx_edges(G,pos=pos,
                           edgelist=edgecoex,
                           edge_color='red',
                           width=edgecoex_size,
                           arrows=True,
                           arrowsize=0.001,
                           alpha=1,
                           node_size=2000,
                           #connectionstyle='Arc3, rad=0.1',
                          ax=ax)
    ##bistable edges##
    nx.draw_networkx_edges(G,pos=pos,
                           edgelist=edgebist,
                           width=edgebist_size,
                           edge_color='blue',
                           arrows=True,
                           arrowsize=0.001,
                           alpha=1,
                           node_size=2000,
                           #connectionstyle='Arc3, rad=0.1',
                          ax=ax)
    ##exclusion edges##
    for edg,edg_size in zip(edgeexclui+edgeexcluj,edgeexclui_size+edgeexcluj_size):
        nx.draw_networkx_edges(G,pos=pos,
                               width=[edg_size*0.5],
                               edgelist=[edg],
                               edge_color='grey',
                               #style='dotted',
                               arrowsize=10*edg_size,
                               arrowstyle='-|>',
                               alpha=0.4,
                               node_size=2000,
                               #connectionstyle="angle3,angleA=-45,angleB=45",
                               ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])

