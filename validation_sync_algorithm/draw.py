import matplotlib.pyplot as plt

from pylab import genfromtxt

import seaborn as sns
sns.set()
#sns.set_context('talk')
sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 1.5})

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

# 0
plot_data.append(genfromtxt("N=100000_fracsize_vs_t.txt"))
labels.append(r'$N=100000$')
colors.append(tableau20[8])

# 1
plot_data.append(genfromtxt("N=10000_fracsize_vs_t.txt"))
labels.append(r'$N=10000$') 
colors.append(tableau20[0])

# 2
plot_data.append(genfromtxt("N=1000_fracsize_vs_t.txt"))
labels.append(r'$N=1000$')
colors.append(tableau20[2])

# 3
plot_data.append(genfromtxt("N=100_fracsize_vs_t.txt"))
labels.append(r'$N=100$')
colors.append(tableau20[4])

# 4
plot_data.append(genfromtxt("N=50_fracsize_vs_t.txt"))
labels.append(r'$N=50$')
colors.append(tableau20[6])

# 5
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.02.txt"))
labels.append(r'ER p=0.02')
colors.append(tableau20[0])

# 6
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.03.txt"))
labels.append(r'ER p=0.03')
colors.append(tableau20[2])

# 7
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.04.txt"))
labels.append(r'ER p=0.04')
colors.append(tableau20[4])

# 9
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.06.txt"))
labels.append(r'ER p=0.06')
colors.append(tableau20[8])

# 10
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.1.txt"))
labels.append(r'ER p=0.1')
colors.append(tableau20[10])

# 11
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.2.txt"))
labels.append(r'ER p=0.2')
colors.append(tableau20[12])

# 12
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.4.txt"))
labels.append(r'ER p=0.4')
colors.append(tableau20[14])

# 13
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.8.txt"))
labels.append(r'ER p=0.8')
colors.append(tableau20[16])

# 14
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=1.txt"))
labels.append(r'ER p=1')
colors.append(tableau20[18])

# 14
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.050.txt"))
labels.append(r'ER p=0.05')
colors.append(tableau20[6])

# 15
plot_data.append(genfromtxt("SF_sigmaVSr_ER_p=0.05.txt"))
labels.append(r'ER p=0.05')
colors.append(tableau20[6])

fig = plt.figure()

nr_plots = len(plot_data)
for i in range(5,13):
    plt.errorbar(100*plot_data[i][:,0], plot_data[i][:,1], plot_data[i][:,2],
                label = labels[i], color = colors[i]
            )    
i=14
plt.errorbar(500*plot_data[i][:,0], plot_data[i][:,1], plot_data[i][:,2],
                label = labels[i], color = colors[i]
            )     

plt.legend(loc="lower right").set_draggable(True)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Formating:
#plt.xlim(0,1.1)
plt.ylim(0,1.1)

#plt.xlabel(r'Edge density ($t$)')
#plt.ylabel(r'$C_{max}/N$')

plt.xlabel(r'Connectivity: ${K/N}$')
plt.ylabel(r'Phase coherence: $r$')

plt.gcf().subplots_adjust(bottom=0.15)
plt.show()