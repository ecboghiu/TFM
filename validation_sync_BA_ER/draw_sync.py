# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:08:11 2019

@author: embog
"""

import matplotlib.pyplot as plt

from pylab import genfromtxt

'''
import seaborn as sns
sns.set()
#sns.set_context('talk')
sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 1.5})
'''

from pylab import rcParams
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Palatino')

golden_ration = (1 + 5 ** 0.5) / 2
one_column_figure_size = 1.7
rcParams['figure.figsize'] = (2*one_column_figure_size * golden_ration*0.9, 2*one_column_figure_size*1.1)
#rcParams['axes.linewidth'] = 0.25
#rcParams['xtick.major.width'] = 0.25
#rcParams['ytick.major.width'] = 0.25


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

plot_data, labels, colors = [], [], []

'''
plot_data.append(genfromtxt("sync_norm_sigmaVSr_ER_N=1000_p=0.006.txt", skip_header=3))
labels.append(r'ER $p$=0.006 $N$=1000')
colors.append(tableau20[2])

plot_data.append(genfromtxt("sync_sigmaVSr_SF_N=1000_gamma=2.27838.txt", skip_header=3))
labels.append(r'SF $\gamma$=2.27838 $N$=1000')
colors.append(tableau20[4])
plot_data.append(genfromtxt("sync_sigmaVSr_SF_N=1000_gamma=1.69285.txt", skip_header=3))
labels.append(r'SF $\gamma$=1.69285 $N$=1000')
colors.append(tableau20[6])
'''
filename = "sync_sigmaVSr_file_N=1000_ER_0.006.txt"
plot_data.append(genfromtxt(filename, skip_header=2))
labels.append(filename)
colors.append(tableau20[0])

filename = "sync_sigmaVSr_file_N=1000_BA_3.txt"
plot_data.append(genfromtxt(filename, skip_header=3))
labels.append(filename)
colors.append(tableau20[6])






fig = plt.figure()
labels[0]=r'Erdős–Rényi (Aleatoria)'
labels[1]=r'Barabási–Albert (Libre de escala)'


i=0
plt.plot((plot_data[i][:,0]), plot_data[i][:,1], #, plot_data[i][:,2],
                label = labels[i],color = colors[i], linewidth=2.0,
            ) 

'''
i=1
plt.plot((plot_data[i][:,0]), plot_data[i][:,1],#, plot_data[i][:,2],
                label = labels[i],color = colors[i], linewidth=2.0,
            ) 
'''

'''
nr_plots = len(plot_data)
for i in range(0,nr_plots):
    plt.errorbar(plot_data[i][:,0], plot_data[i][:,1], plot_data[i][:,2],
                label = labels[i], color = colors[i]
            )
'''

plt.legend(loc="upper left").set_draggable(True)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Formating:
plt.xlim(0,0.2)
plt.ylim(0,1)

#plt.xlabel(r'Edge density ($t$)')
#plt.ylabel(r'$C_{max}/N$')

plt.xlabel(r'Acoplo $K/N$')
plt.ylabel(r'Coherencia $r$')

plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('sync_BA_ER1.png',dpi=1200)
plt.show()