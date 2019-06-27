# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 09:36:05 2019

@author: embog
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.integrate as inte
from matplotlib.colors import LinearSegmentedColormap

'''
from pylab import rcParams
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Palatino')
'''


#para que en jupytter se vean mejor las graficas
#get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")

# folllowing lines taken form: https://github.com/cbnfreitas/kuramoto_model_integrate_and_plot
def chop_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#This cmap is applied in to node coloring 
#cmap_aux = LinearSegmentedColormap.from_list("", ["cyan","#580058"])
#see https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib for colors
cmap_aux = LinearSegmentedColormap.from_list("", ["red","yellow"])

#In gray scale, cyan becomes almost white, this is why we chop the begining of the color map
cyan_purple_cmap = chop_colormap(cmap_aux, 0.00, 1)
ccmap_aux = LinearSegmentedColormap.from_list("", ["darkred","red","yellow","red","darkred"])
'''
viridis = plt.cm.get_cmap('viridis', 12)
ccmap_aux = LinearSegmentedColormap.from_list("", [viridis(.01),viridis(.1),viridis(.2),
                                                   viridis(.3),viridis(.4),viridis(.5),
                                                   viridis(.6),viridis(.7),viridis(.8),viridis(0.9),viridis(.99),
                                                  viridis(.9),viridis(.8),viridis(.7),
                                                   viridis(.6),viridis(.5),viridis(.4),
                                                  viridis(.3),viridis(.2),viridis(.1),
                                                  viridis(.01),])
'''
psi_purple_cmap = chop_colormap(ccmap_aux, 0.00, 1)

def frequency_to_color(w):
    colors_cyan_purple = cyan_purple_cmap(np.linspace(0, 1, 1001))
    w_min = min(w)
    w_max = max(w)
    return [colors_cyan_purple[int(1000*(w[i] - w_min)/(w_max - w_min))] for i in range(len(w))]

def angle_to_color(w):
    colors_cyan_purple = psi_purple_cmap(np.linspace(0, 1, 1001))
    w_min = min(w)
    w_max = max(w)
    return [colors_cyan_purple[int(1000*(w[i] - w_min)/(w_max - w_min))] for i in range(len(w))]


# me gustan más estos colores para los plots
def colores_tableau():
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
        red, green, blue = tableau20[i]    
        tableau20[i] = (red / 255., green / 255., blue / 255.) 
    return tableau20

tableau20 = colores_tableau()

filename_in = "../data/FG_N=500_m=1_4_a=5_sig=0.08_EDGELIST_1.txt"
datos = (np.genfromtxt(filename_in,skip_header=0))
xdat=datos[:,1]
ydat=datos[:,2]
e1 = [int(xdat[i]) for i in range(0,len(xdat))]
e2 = [int(ydat[i]) for i in range(0,len(ydat))]
my_edge_list = [[e1[i],e2[i],datos[i,0]] for i in range(0,len(e1))]
my_edge_list = np.array(my_edge_list, dtype=np.int)
I1 = np.max(my_edge_list[:,0])
I2 = np.max(my_edge_list[:,1])
NODE_NR = np.max([I1,I2])+1
NODE_NR


# In[19]:


data = np.genfromtxt("../pattern.txt")
data = np.transpose(data)
list_dat = list(data)
t = list_dat.pop(0)
delta_t_trash = list_dat.pop(0)
asd = np.array(list_dat)
t = [int(t[i]) for i in range(len(t))]

nr_e_ini = np.min(t)
nr_e_fin = np.max(t)
nr_e_ini
    


# In[64]:


g_aux = nx.Graph()
for i in range(0,NODE_NR):
    g_aux.add_node(i)
for i in range(0,len(my_edge_list)):
    g_aux.add_edge(my_edge_list[i][0],my_edge_list[i][1])

g_aux.add_edge(0,432)
    
#connctd = sorted(nx.connected_components(g_aux), key=len, reverse=True)
#connctd = [np.array((list(i)),dtype=int) for i in connctd]
    
fig = plt.figure(figsize=(11,9))
POS_FINAL = nx.spring_layout(g_aux, iterations=40)
nx.draw(g_aux, pos=POS_FINAL, node_size=40,
        edge_color = 'black',
        with_labels=False, font_size=10, font_color='white')


# Ahora vamos a animar.

# In[65]:


import matplotlib.animation as animation
from IPython.core.display import display, HTML
display(HTML('<h1> Animation attempt. </h1>'))


import matplotlib.animation as animation

#%pylab qt
 

plt.rcParams['animation.ffmpeg_path'] = 'c://ffmpeg-4.1.1-win64-static//bin//ffmpeg.exe'
fig = plt.figure(figsize=(11,9))
#param, = plt.plot(t,y_t[1])
ang_2pi = np.linspace(0,2*np.pi,num=1000)
#colors_w = frequency_to_color(w)


def update_argand(i):
    fig.clf()
    
    
    for node in range(0,NODE_NR):
        plt.scatter(np.cos(asd[node][i]),np.sin(asd[node][i]),facecolor=colors_w[node], # color='b',
                marker='o', linewidth=0.5, zorder=3, edgecolor='gray', s = 100)
        
    plt.scatter(r[i]*np.cos(psi[i]),r[i]*np.sin(psi[i]), color='b',
                linewidth = 3, zorder=2)
    plt.plot(r[:i]*np.cos(psi[:i]),r[:i]*np.sin(psi[:i]), color='b', linewidth = 1.0, zorder=1)
    
    plt.plot(np.cos(ang_2pi),np.sin(ang_2pi), color='gray')
    
    plt.xlim([-1.1,1.1])
    plt.ylim([-1.1,1.1])
    
    sm = plt.cm.ScalarMappable(cmap=cyan_purple_cmap, norm=plt.Normalize(vmin=min(w), vmax=max(w)))
    sm._A = []
    cb1 = plt.colorbar(sm)
    cb1.set_label(r'$\omega_i$')
    
    return fig,

g = nx.Graph()
for i in range(0,NODE_NR):
    g.add_node(i)

i = 0
while (i<nr_e_ini):
    g.add_edge(my_edge_list[i][0],my_edge_list[i][1])
    i = i + 1
    
w = asd[:,0]
colors_s = angle_to_color(w)
nx.draw(g, pos=POS_FINAL, node_size=40, node_color = colors_s,
        edge_color = 'black',
        with_labels=False, font_size=10, font_color='white')

#my_edge_list = list(g.edges())

asd = np.mod(asd, 2*np.pi)
y_t=asd

number_of_edges = g.number_of_edges()

    
edge_to_print = nr_e_ini
edge_n1_toprint = my_edge_list[edge_to_print][0]
edge_n2_toprint = my_edge_list[edge_to_print][1]

if edge_to_print != t[0]:
    print("warning 2")
    
def update_graph(i):
    fig.clf()
    
    w = asd[:,i]
    colors_s = angle_to_color(w)
    nx.draw(g, pos=POS_FINAL,
            node_size=40, node_color = colors_s,
        edge_color = 'black',
        with_labels=False, font_size=10, font_color='white')
    
    ax = plt.gca()
    ax.collections[0].set_edgecolor("black") #Drawing a black border around nodes.
    ax.collections[0].set_linewidth(0.5)
    ax.collections[1].set_linewidth(0.5)  # Change edges linewidth
    sm = plt.cm.ScalarMappable(cmap=psi_purple_cmap, norm=plt.Normalize(vmin=0, vmax=2*np.pi))
    sm._A = []
    cb1 = plt.colorbar(sm)
    cb1.set_label(r'$\psi$')
    
    edge_to_print = t[i]
    edge_n1_toprint = my_edge_list[edge_to_print][0]
    edge_n2_toprint = my_edge_list[edge_to_print][1]
        
    g.add_edge(edge_n1_toprint, edge_n2_toprint)
    H = g.edge_subgraph([(edge_n1_toprint, edge_n2_toprint)])
    nx.draw_networkx_edges(H, pos=POS_FINAL,
                       with_labels=False,
                       edge_color='r',
                       width=7.0,
                       alpha=1.0,
                       font_size=10,
                      )
    #nx.draw(g, pos=POS_FINAL, node_size=20, #node_color = colors_w,
    #    edge_color = 'black',
    #    with_labels=False, font_size=10, font_color='white')
    
    
    '''
    for node in range(0,NODE_NR):
        plt.scatter(np.cos(y_t[node][i]),np.sin(y_t[node][i]),facecolor=colors_w[node], # color='b',
                marker='o', linewidth=0.5, zorder=3, edgecolor='gray', s = 100)
        
    plt.scatter(r[i]*np.cos(psi[i]),r[i]*np.sin(psi[i]), color='b',
                linewidth = 3, zorder=2)
    plt.plot(r[:i]*np.cos(psi[:i]),r[:i]*np.sin(psi[:i]), color='b', linewidth = 1.0, zorder=1)
    
    plt.plot(np.cos(ang_2pi),np.sin(ang_2pi), color='gray')
    
    plt.xlim([-1.1,1.1])
    plt.ylim([-1.1,1.1])
    
    sm = plt.cm.ScalarMappable(cmap=cyan_purple_cmap, norm=plt.Normalize(vmin=min(w), vmax=max(w)))
    sm._A = []
    cb1 = plt.colorbar(sm)
    cb1.set_label(r'$\omega_i$')
    '''
    plt.title("$\ell$=",t[i])
    return fig,

# this if for the frame by frame kuramoto plot
'''
i=1600
fig.clf()
for node in range(0,NODE_NR):
    plt.scatter(np.cos(asd[node][i]),np.sin(asd[node][i]),facecolor=colors_w[node], # color='b',
            marker='o', linewidth=0.5, zorder=3, edgecolor='gray', s = 200)

plt.scatter(r[i]*np.cos(psi[i]),r[i]*np.sin(psi[i]), color='b',
            linewidth = 3, zorder=4)
#plt.plot(r[:i]*np.cos(psi[:i]),r[:i]*np.sin(psi[:i]), color='b', linewidth = 1.0, zorder=1)

plt.plot(np.cos(ang_2pi),np.sin(ang_2pi), color='gray',zorder=1)

plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])

sm = plt.cm.ScalarMappable(cmap=cyan_purple_cmap, norm=plt.Normalize(vmin=min(w), vmax=max(w)))
sm._A = []
cb1 = plt.colorbar(sm)
cb1.set_label(r'$\omega_i$')
plt.savefig("instant_BI_2.pdf")
'''



frame = [i for i in range(0,len(t))]
ani = animation.FuncAnimation(fig, update_graph, frames=frame,
                              interval=100, blit=True)

import datetime

currentDT = datetime.datetime.now()
print (str(currentDT))
#ani = animation.FuncAnimation(fig, update_argand, frames=10,interval=10, blit=True)

#ani.save('basic_animation.mp4', fps=60, extra_args=['-vcodec', 'libx264'])
#ani.save('test.gif', writer='imagemagick')

#ani.save("movie.mp4")
HTML(ani.to_html5_video())

currentDT = datetime.datetime.now()
print (str(currentDT))

#plt.show()


# In[ ]:





# In[66]:


HTML(ani.to_html5_video())


# In[ ]:


#ani.save('network_chula2.mp4', fps=60, extra_args=['-vcodec', 'libx264'])


# In[ ]:

