# -*- coding: utf-8 -*-
"""
Created on Tue May 28 15:47:38 2019

@author: embog
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import networkx as nx
import collections

datos = (np.genfromtxt("../data/FG_N=500_m=1_a=10_sig=0.08_EDGELIST.txt",
                       skip_header=0))

xdat=datos[:,1]
ydat=datos[:,2]


NODE_NR = 500
G=nx.Graph()
for i in range(0,NODE_NR):
    G.add_node(i)

options = {
    'node_color': 'black',
    'node_size': 50,
    'line_color': 'grey',
    'linewidths': 0,
    'width': 0.1,
}

hist_surf = []
for x, y in zip(xdat,ydat):
    G.add_edge(x,y)
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    hist_surf.append([deg, cnt])
    
    
aux = np.linspace(0,len(hist_surf)/NODE_NR,len(hist_surf),dtype=np.float32)
aux2 = np.ones(len(hist_surf))
fig = plt.figure()
#ax = Axes3D(fig)
file = open("hist3D.txt","w")
for i in range(0,len(hist_surf)):
    for j in range(0,len(hist_surf[i][0])):   
        file.write(str(aux[i])+" "+str(hist_surf[i][0][j])+" "+str(hist_surf[i][1][j])+"\n")
file.close()       
    

    
    