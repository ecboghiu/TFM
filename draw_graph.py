#print("Importing plot.py!")

#from __future__ import division # so that 1/2 does not give the loathed 0
#from IPython.display import display, Math, Latex # beautiful output
import matplotlib.pyplot as plt
import collections
import numpy as np
import seaborn as sns
sns.set()
#sns.set_context('talk')
sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 1.5})
import networkx as nx



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

NODE_NR = 1000
GAMMA = 1.69285028
K_MIN = 2
K_MAX = NODE_NR


def promedio_dif (y):
    k_avg = 6
    
    K_MIN = 2
    K_MAX = int(NODE_NR**0.5)
    
    suma = 0
    for k in range(K_MIN,K_MAX):
        suma = suma + k**(-y)
        
    suma2 = 0
    for k in range(K_MIN,K_MAX):
        suma2 = suma2 + k* (k**(-y))/suma
        
    return k_avg - suma2 / suma, suma, suma2

NORMALIZING_FACTOR = 1/(promedio_dif(GAMMA)[1])
print("Norm const python:", NORMALIZING_FACTOR)
print("promedio teorico:", (promedio_dif(GAMMA)[2]))

G = nx.read_adjlist("adj_C.txt", nodetype=int)

options = {
    'node_color': 'black',
    'node_size': 50,
    'line_color': 'grey',
    'linewidths': 0,
    'width': 0.1,
}
#nx.draw(G, **options)
#plt.show()

degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
# print "Degree sequence", degree_sequence
degreeCount = collections.Counter(degree_sequence)
deg, cnt = zip(*degreeCount.items())

average = 0
norm = 0
for i in range(0,len(deg)):
    norm = norm + cnt[i]    
for i in range(0,len(deg)):
    average = average + deg[i] * cnt[i]
average = average/norm
print("average=",average)

fig, ax = plt.subplots()
plt.plot(np.array(deg), np.array(cnt)/norm, color='b')#, width=0.80, color='b')
plt.plot(np.array(deg),
         NORMALIZING_FACTOR * np.array(deg)**(-GAMMA),
   label = 'teor', color = tableau20[1])
    
#plt.yscale('log')
#plt.xscale('log')

plt.title("Degree Histogram")
plt.ylabel("Count")
plt.xlabel("Degree")
ax.set_xticks([d + 0.4 for d in deg])
ax.set_xticklabels(deg)

# draw graph in inset
plt.axes([0.4, 0.4, 0.5, 0.5])
Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]
pos = nx.spring_layout(G)
plt.axis('off')
nx.draw_networkx_nodes(G, pos, node_size=20)
nx.draw_networkx_edges(G, pos, alpha=0.4)

plt.show()

