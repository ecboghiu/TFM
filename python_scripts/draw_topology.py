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

def hist_3d (filename_in, filename_out):
    datos = (np.genfromtxt(filename_in,
                           skip_header=5))
    xdat=datos[:,1]
    ydat=datos[:,2]
    NODE_NR = 2000
    G=nx.Graph()
    for i in range(0,NODE_NR):
        G.add_node(i)
    '''
    options = {
        'node_color': 'black',
        'node_size': 50,
        'line_color': 'grey',
        'linewidths': 0,
        'width': 0.1,
    }
    '''
    hist_surf = []
    for x, y in zip(xdat,ydat):
        G.add_edge(x,y)
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = zip(*degreeCount.items())
        
        np_deg = np.array(deg)
        np_cnt = np.array(cnt)
        norm = np.sum(np_cnt)   
        np_prob = np_cnt/norm
        average = np.sum(np_deg*np_prob)
        var = np.sum(np_deg*np_deg*np_prob)-average
        
        hist_surf.append([deg, cnt, average, var])
        
    fig = plt.figure()
    i=10
    plt.xlabel(r"$\mathrm{log} k$")
    plt.ylabel(r"$\mathrm{log} f$")
    #for j in range(0,len(hist_surf)): 
    #j=(len(hist_surf)-1)
    #plt.plot((hist_surf[j][0][:]),(hist_surf[j][1][:]))
    

    
    aux = np.linspace(0,len(hist_surf)/NODE_NR,len(hist_surf),dtype=np.float32)
    aux2 = np.ones(len(hist_surf))
    
    
    #ax = fig.add_subplot(111, projection='3d')
    #ax = Axes3D(fig)
    '''
    file = open(filename_out,"w")
    for i in range(0,len(hist_surf)):
        #if (i % (len(hist_surf)/3) == 0 ):
        for j in range(0,len(hist_surf[i][0])):   
            if (j % (len(hist_surf[i][0])/(len(hist_surf[i][0])/10)) == 0 ):
                file.write(str(aux[i])+" "+str(np.log(hist_surf[i][0][j]))+" "+str(np.log(hist_surf[i][1][j]))+"\n")
    file.close()  
    '''    
    
    print("End of function"+filename_in+filename_out)
    
alpha="5"
sig = "0.01"
file_in = "../data/FG_N=500_m=1_a="+alpha+"_sig="+sig+"_EDGELIST.txt"
file_out = "../data/hist3D.txt"
hist_3d(file_in,file_out)

'''    
lam=['0.05','0.08','0.1','0.5','1']
alpha=['-20','-10','-5','-1','0','1','5','10','20']
for l in lam:
    for a in alpha:
        file_in="../data/FG_N=2000_m=1_a="+a+"_sig="+l+"_EDGELIST.txt"
        file_out = "../data/hist3D_2000_"+a+"_"+l+".txt"
        hist_3d(file_in,file_out)
'''