#print("Importing plot.py!")

#from __future__ import division # so that 1/2 does not give the loathed 0
#from IPython.display import display, Math, Latex # beautiful output
import matplotlib.pyplot as plt

import numpy as np
import seaborn as sns


#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True

#### PLOT
    
# Sources: http://www.randalolson.com/2014/06/28/how-to-make
#            -beautiful-data-visualizations-in-python-with-matplotlib/
# These are the "Tableau 20" colors as RGB.   
tableau20=[(31, 119, 180),  (174, 199, 232), (255, 127, 14),  (255, 187, 120),    
            (44, 160, 44),   (152, 223, 138), (214, 39, 40),   (255, 152, 150),    
            (148, 103, 189), (197, 176, 213), (140, 86, 75),   (196, 156, 148),    
            (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
            (188, 189, 34),  (219, 219, 141), (23, 190, 207),  (158, 218, 229)]
# Scale the RGB values to the [0, 1] range, which is the format matplotlib
# accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 

def plot_hist():
    plot_data = []
    input_data = open('hist.txt','r')
    for line in input_data.readlines():
        aux_line = (line.strip()).split()
        aux_list = []
        for e in range(0,2):
            if aux_line[0] == '#':
                break
            aux_list.append(float(aux_line[e]))
        plot_data.append(np.array(aux_list))
    input_data.close()

    x, y = [], []
    for i_aux in range(1,len(plot_data)):
        x.append(plot_data[i_aux][0])
        y.append(plot_data[i_aux][1])
    x = np.array(x)
    y = np.array(y)

    labels = [r'$c_{q}$', r'$c_{p}$']
    plt.figure(figsize=(6, 7)) 
    plt.xlim(0,1.2)
    plt.style.use('ggplot')
    #plt.style.use('seaborn-talk')
    plt.xlabel(r'Degree')
    plt.ylabel(r'Freq')
    plt.ylim(0,1.2)
    for el in range(0,1):    
        plt.plot(x, y,  color=tableau20[el], label = labels[el])
    plt.title(r'Histogram')
    plt.legend(loc='upper right')
    plt.show()

def plot_sigma_vs_r():
    plot_data = []
    input_data = open('SF_sigmaVSr.txt','r')
    for line in input_data.readlines():
        aux_line = (line.strip()).split()
        aux_list = []
        for e in range(0,3):
            if aux_line[0] == '#':
                break
            aux_list.append(float(aux_line[e]))
        plot_data.append(np.array(aux_list))
    input_data.close()

    x, y, yerr= [], [], []
    for i_aux in range(1,len(plot_data)):
        x.append(plot_data[i_aux][0])
        y.append(plot_data[i_aux][1])
        yerr.append(plot_data[i_aux][2])
    x = np.array(x)
    y = np.array(y)
    yerr = np.array(yerr)

    labels = [r'$c_{q}$', r'$c_{p}$']
    plt.figure(figsize=(8, 7)) 
    plt.xlim(0,1.1)
    plt.style.use('ggplot')
    plt.style.use('seaborn-talk')
    plt.xlabel(r'$\sigma$')
    plt.ylabel(r'$r$')
    plt.ylim(0,1.1)
    for el in range(0,1):    
        plt.plot(x, y,  color=tableau20[el], label = labels[el])
        #ax = plt.gca()
        #ax.errorbar(x, y, yerr=yerr, fmt='o',
        #                color=tableau20[el], label = labels[el])
    plt.title(r'Phase transition')
    plt.legend(loc='upper right')
    plt.show()

def plot_perc():
    plt.figure()
    sns.set()
    sns.set_context('talk')
    sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 2})

    plot_data = []
    input_data = open('frac_size_vs_t.txt','r')
    for line in input_data.readlines():
        aux_line = (line.strip()).split()
        aux_list = []
        for e in range(0,len(aux_line)):
            aux_list.append(float(aux_line[e]))
        plot_data.append(np.array(aux_list))
    input_data.close()

    
    x, y= [], []
    for i_aux in range(1,len(plot_data)):
        x.append(plot_data[i_aux][0])
        y.append(plot_data[i_aux][1])
    x = np.array(x)
    y = np.array(y)

    labels = [r'$C_{max}/N$']
    #plt.figure(figsize=(7, 5)) 
    plt.xlim(0,1.1)

    plt.xlabel(r'Edge density ($t$)')
    plt.ylabel(r'$C_{max}/N$')
    plt.ylim(0,1.1)
    for el in range(0,1):    
        plt.plot(x, y,  color=tableau20[0], label = labels[el])

    plt.title(r'Network with product rule.')
    plt.legend(loc='upper left')

    #plt.tight_layout() #to make room for xlabel
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.show()

def plot_graph(filename):
    import networkx as nx
    plt.figure()
    G = nx.read_adjlist(filename, nodetype=int)
    #plt.figure()
    options = {
        'node_color': 'black',
        'node_size': 50,
        'line_color': 'grey',
        'linewidths': 0,
        'width': 0.1,
    }
    nx.draw(G, **options)
    '''
    # taken from https://networkx.github.io/documentation/
    # stable/auto_examples/drawing/plot_giant_component.htm
    # l#sphx-glr-auto-examples-drawing-plot-giant-component-py
    
    #region = 220
    #plt.subplots_adjust(left=0, right=1, bottom=0, top=0.95,
    #                                wspace=0.01, hspace=0.01)
    
    #import math
    #try:
    #    import pygraphviz
    #    from networkx.drawing.nx_agraph import graphviz_layout
    #    layout = graphviz_layout
    #except ImportError:
    #    try:
    #        import pydot
    #        from networkx.drawing.nx_pydot import graphviz_layout
    #        layout = graphviz_layout
    #    except ImportError:
    #        print("PyGraphviz and pydot not found;\n"
    #            "drawing with spring layout;\n"
    #            "will be slow.")
    #        layout = nx.spring_layout
    
    layout = nx.spring_layout
    pos = layout(G)
    #region += 1
    #plt.subplot(region)
    nx.draw(G, pos,
            with_labels=False,
            node_size=10
            )
    # identify largest connected component
    Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
    G0 = Gcc[0]
    nx.draw_networkx_edges(G0, pos,
                            with_labels=False,
                            edge_color='r',
                            width=6.0
                            )
    # show other connected components
    for Gi in Gcc[1:]:
        if len(Gi) > 1:
            nx.draw_networkx_edges(Gi, pos,
                                    with_labels=False,
                                    edge_color='r',
                                    alpha=0.3,
                                    width=5.0
                                    )
    
    '''
    plt.show()



#plot_graph("adj_C.txt")
#plt.figure()
import networkx as nx
G = nx.read_adjlist("adj_C.txt", nodetype=int)
options = {
    'node_color': 'black',
    'node_size': 50,
    'line_color': 'grey',
    'linewidths': 0,
    'width': 0.1,
}
nx.draw(G, **options)
plt.show()
#Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
#print("largest connected:", Gcc[0].size())

#plot_perc()

#plot_hist()
#plot_sigma_vs_r()

#plot_graph_visually('adj_C.txt')
'''
from matplotlib import pyplot
from pylab import genfromtxt
mat0 = genfromtxt("N=10000_fracsize_vs_t.txt")
#mat1 = genfromtxt("data1.txt");
pyplot.plot(mat0[:,0], mat0[:,1], label = "data0")
#pyplot.plot(mat1[:,0], mat1[:,1], label = "data1")
pyplot.legend()
pyplot.xlim(0,1.1)
pyplot.xlabel(r'Edge density ($t$)')
pyplot.ylabel(r'$C_{max}/N$')
pyplot.ylim(0,1.1)
pyplot.show()
'''