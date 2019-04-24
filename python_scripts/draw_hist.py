import matplotlib.pyplot as plt

from pylab import genfromtxt
import numpy as np

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
plot_data.append(genfromtxt("hist.txt", skip_header=1))
labels.append(r'ER')
colors.append(tableau20[0])

fig = plt.figure()

NODE_NR = len(plot_data[0][:,0])
GAMMA = 2.27837758
K_MIN = 2
K_MAX = NODE_NR


def promedio_dif (y):
    k_avg = 6
    
    K_MIN = 2
    K_MAX = NODE_NR
    
    suma = 0
    for k in range(K_MIN,K_MAX):
        suma = suma + k**(-y)
        
    suma2 = 0
    for k in range(K_MIN,K_MAX):
        suma2 = suma2 + k* k**(-y)
        
    return k_avg - suma2 / suma, suma, suma2, suma2 / suma

NORMALIZING_FACTOR = 1/(promedio_dif(GAMMA)[1])
print("Norm const python:", NORMALIZING_FACTOR)
print("promedio teorico:", (promedio_dif(GAMMA)[3]))

nr_plots = len(plot_data)
for i in [0]:#range(5,13):
    plt.plot(plot_data[i][:,0], plot_data[i][:,1],
            label = labels[i], color = colors[0]
            )
    #plt.plot(   np.array(plot_data[i][:,0]),
    #        NORMALIZING_FACTOR * np.array(plot_data[i][:,0])**(-GAMMA) ,
    #            label = 'Thero', color = tableau20[1]
    #        )


#plt.legend(loc="lower right").set_draggable(True)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Formating:
#plt.xlim(0,1.1)
#plt.ylim(0,1.1)


plt.xlabel(r'$K$')
plt.ylabel(r'Relative frequency')

plt.title('Degree histogram')

plt.gcf().subplots_adjust(bottom=0.15)
plt.show()