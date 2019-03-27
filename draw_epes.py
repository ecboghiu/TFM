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


filename = r"trampa_EPES_sigmaVSr_file_N=10000_weff.txt"
plot_data.append(genfromtxt(filename,  skip_header=2))
labels.append(['r_'+filename,'p_'+filename])
colors.append([tableau20[4],tableau20[2]])

fig = plt.figure()

for i in range(0,len(plot_data)):
    plt.errorbar((plot_data[i][:,0]), plot_data[i][:,2], plot_data[i][:,3],
                label = labels[i][0], color = colors[i][0])
    plt.plot((    plot_data[i][:,0]), plot_data[i][:,1],
                label = labels[i][1], color = colors[i][1]) 

plt.legend().set_draggable('True')

#plt.yscale('log')
# Formating:
#plt.xlim(0,1.1)
plt.xlabel(r'Edge density ($t$)')
plt.ylabel(r'$C_{max}/N$')
plt.ylim(0,1.1)
plt.title(r'Synchronization and maximum component fraction')

plt.gcf().subplots_adjust(bottom=0.15)
plt.show()