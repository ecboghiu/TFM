import matplotlib.pyplot as plt

from pylab import genfromtxt
import numpy as np


from pylab import rcParams
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Palatino')


golden_ration = (1 + 5 ** 0.5) / 2
one_column_figure_size = 1.7
rcParams['figure.figsize'] = (2*one_column_figure_size * golden_ration, 2*one_column_figure_size)
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

r'''
# 0
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.02.txt",skip_header=2))
labels.append(r'ER p=0.02')
colors.append(tableau20[1])

# 1
plot_data.append(genfromtxt("SF_sigmaVSr_SF_gamma=2.21412.txt", skip_header=2))
labels.append(r'SF $\gamma$=2.21412')
colors.append(tableau20[3])

plot_data.append(genfromtxt("coh_0.7_300.txt", skip_header=2))
labels.append(r'$\sigma$=0.7')
colors.append(tableau20[8])


plot_data.append(genfromtxt("coh_0.8_300.txt", skip_header=2))
labels.append(r'$\sigma$=0.8')
colors.append(tableau20[10])
'''
'''
plot_data.append(genfromtxt("coh_1e-005.txt", skip_header=2))
labels.append(r'$h$=0.00001')
colors.append(tableau20[10])
plot_data.append(genfromtxt("coh_0.0001.txt", skip_header=2))
labels.append(r'$h$=0.0001')
colors.append(tableau20[0])
plot_data.append(genfromtxt("coh_0.001.txt", skip_header=2))
labels.append(r'$h$=0.001')
colors.append(tableau20[2])
plot_data.append(genfromtxt("coh_0.01.txt", skip_header=2))
labels.append(r'$h$=0.01')
colors.append(tableau20[4])
plot_data.append(genfromtxt("coh_0.1.txt", skip_header=2))
labels.append(r'$h$=0.1')
colors.append(tableau20[6])
plot_data.append(genfromtxt("coh_0.5.txt", skip_header=2))
labels.append(r'$h$=0.5')
colors.append(tableau20[8])
'''
plot_data.append(genfromtxt("../data/FG_N=500_m=1_a=0_sig=1.5_h_high.txt", skip_header=5))
labels.append(r'$h$=0.01')
colors.append(tableau20[0])

plot_data.append(genfromtxt("../data/FG_N=500_m=1_a=0_sig=1.5_fantastico.txt", skip_header=5))
labels.append(r'$h$=0.1')
colors.append(tableau20[2])

plot_data.append(genfromtxt("../data/FG_N=500_m=1_a=0_sig=1.5_hvlow.txt", skip_header=5))
labels.append(r'$h$=0.5')
colors.append(tableau20[4])

plot_data.append(genfromtxt("../data/FG_N=500_m=1_a=0_sig=1.5_nlow.txt", skip_header=5))
labels.append(r'$h$=0.1 low')
colors.append(tableau20[4])

fig = plt.figure()

nr_plots = len(plot_data)
'''
for i in [0]:
    plt.errorbar(plot_data[i][:,0], plot_data[i][:,1], plot_data[i][:,2],
                label = labels[i], color = colors[i]
            )   
'''
for i in range(0,len(plot_data)):
    plt.plot(plot_data[i][:,0],
             (plot_data[i][:,2]),
                label = labels[i], color = colors[i])
'''
plt.plot(plot_data[i][:,0], plot_data[i][:,1], 
                label = labels[i], color = tableau20[0])  
plt.plot(plot_data[i][:,0], plot_data[i][:,2],
                label = 'teor', color = tableau20[2])  
'''


  

#plt.yscale('log')
plt.legend(loc="upper left").set_draggable(True)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Formating:
#plt.xlim(0,1.1)
#plt.ylim(0,1.1)

#plt.xlabel(r'Edge density ($t$)')
#plt.ylabel(r'$C_{max}/N$')

plt.xlabel(r'$\ell$')
plt.ylabel(r'r')
plt.title(r'Variación en $r$ debido a cambios en $h$')

plt.gcf().subplots_adjust(bottom=0.15)
plt.show()