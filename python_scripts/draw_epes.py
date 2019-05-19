import matplotlib.pyplot as plt

import numpy as np
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
data_labels_colors = [plot_data, labels, colors]

def take_from_txt (out, filename, color):
    out[0].append(genfromtxt(filename,  skip_header=5))
    out[1].append(['r_'+filename,'p_'+filename])
    color_code = color
    out[2].append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 
   
def take_from_txt_Nas (out, N, a, sig, color):
    filename = r"../data/FG_N="+str(N)+"_m=1_a="+str(a)+"_sig="+str(sig)+".txt"
    out[0].append(genfromtxt(filename,  skip_header=5))
    out[1].append(['r_'+filename,'p_'+filename])
    color_code = color
    out[2].append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 
  

# Sigma slices
def plot_sig_slice(SIG):
    #take_from_txt_Nas(data_labels_colors, 200, -3, SIG, 14)
    take_from_txt_Nas(data_labels_colors, 200, -1, SIG, 0)
    take_from_txt_Nas(data_labels_colors, 200, -0.5, SIG, 2)
    take_from_txt_Nas(data_labels_colors, 200, 0, SIG, 4)
    take_from_txt_Nas(data_labels_colors, 200, 0.5, SIG, 6)
    take_from_txt_Nas(data_labels_colors, 200, 1, SIG, 8)
    take_from_txt_Nas(data_labels_colors, 200, 3, SIG, 10)
    take_from_txt_Nas(data_labels_colors, 200, 5, SIG, 12)
    
# Alpha slices
def plot_alph_slice (ALPHA):
    take_from_txt_Nas(data_labels_colors, 200, ALPHA, 0.5, 0)
    take_from_txt_Nas(data_labels_colors, 200, ALPHA, 1, 2)
    take_from_txt_Nas(data_labels_colors, 200, ALPHA, 1.5, 4)
    take_from_txt_Nas(data_labels_colors, 200, ALPHA, 3, 6)

#plot_sig_slice(0.1)

#plot_alph_slice(0)

#ALPHA=-1
#take_from_txt_Nas(data_labels_colors, 200, 0, 1.5, 0)
'''
take_from_txt_Nas(data_labels_colors, 1000, -1, "1.50", 0)
take_from_txt_Nas(data_labels_colors, 1000, 0, "1.50", 2)
take_from_txt_Nas(data_labels_colors, 1000, 1, "1.50", 4)
'''


'''
take_from_txt_Nas(data_labels_colors, 1000, -0.8, 1.5, 0)
take_from_txt_Nas(data_labels_colors, 1000, 0, 1.5, 2)
take_from_txt_Nas(data_labels_colors, 1000, +0.8, 1.5, 4)
'''


take_from_txt_Nas(data_labels_colors, 100, 0.5, 3, 0)


'''
filename = r"../data/FG_N=200_m=1_a=-0.5_sig=0.1.txt"
plot_data.append(genfromtxt(filename,  skip_header=5))
labels.append(['r_'+filename,'p_'+filename])
color_code = 0
colors.append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 

filename = r"../data/FG_N=200_m=1_a=-0.5_sig=1.txt"
plot_data.append(genfromtxt(filename,  skip_header=5))
labels.append(['r_'+filename,'p_'+filename])
color_code = 2
colors.append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 

filename = r"../data/FG_N=200_m=1_a=-0.5_sig=1.5.txt"
plot_data.append(genfromtxt(filename,  skip_header=5))
labels.append(['r_'+filename,'p_'+filename])
color_code = 4
colors.append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 
'''

'''
filename = r"../data/FG_N=200_m=1_a=-0.5_sig=0.1.txt"
plot_data.append(genfromtxt(filename,  skip_header=5))
labels.append(['r_'+filename,'p_'+filename])
color_code = 2
colors.append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 


filename = r"../data/FG_N=200_m=1_a=0.5_sig=0.1.txt"
plot_data.append(genfromtxt(filename,  skip_header=5))
labels.append(['r_'+filename,'p_'+filename])
color_code = 4
colors.append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]]) 
'''


fig = plt.figure()

for i in range(0,len(plot_data)):
    plt.errorbar((plot_data[i][:,0]), plot_data[i][:,2], np.array(plot_data[i][:,3])/np.sqrt(3000),
                label = labels[i][0], color = colors[i][0])
    plt.plot((    plot_data[i][:,0]), plot_data[i][:,1],
                label = labels[i][1], color = colors[i][1])
    #plt.plot((    plot_data[i][:,0]), plot_data[i][:,4],
    #            label = labels[i][1], color = colors[i][2])

plt.legend().set_draggable('True')

#plt.yscale('log')
# Formating:
#plt.xlim(0,1.1)
plt.xlabel(r'Edge density ($t$)')
plt.ylabel(r'$C_{max}/N$,r')
plt.ylim(0,1.1)
plt.title(r'Synchronization and maximum component fraction')

plt.gcf().subplots_adjust(bottom=0.15)
plt.show()