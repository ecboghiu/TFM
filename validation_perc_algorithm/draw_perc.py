import matplotlib.pyplot as plt

from pylab import genfromtxt


import seaborn as sns
#sns.set()
#sns.set_context('talk')
#sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 1.5})


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

'''
plot_data.append(genfromtxt("N=100000_fracsize_vs_t_old.txt"))
labels.append(r'$N=100000$')
colors.append(tableau20[8])

plot_data.append(genfromtxt("N=10000_fracsize_vs_t_old.txt"))
labels.append(r'$N=10000$') 
colors.append(tableau20[0])

plot_data.append(genfromtxt("N=1000_fracsize_vs_t_old.txt"))
labels.append(r'$N=1000$')
colors.append(tableau20[2])

plot_data.append(genfromtxt("N=100_fracsize_vs_t_old.txt"))
labels.append(r'$N=100$')
colors.append(tableau20[4])

plot_data.append(genfromtxt("N=50_fracsize_vs_t_old.txt"))
labels.append(r'$N=50$')
colors.append(tableau20[6])
'''

#plot_data.append(genfromtxt("N=500000_fracsize_vs_t.txt"))
#labels.append(r'$N=500000$ Null')
#colors.append(tableau20[6])

plot_data.append(genfromtxt("N=100000_fracsize_vs_t.txt"))
labels.append(r'$N=10000$ Aleatorio')
colors.append(tableau20[0])

plot_data.append(genfromtxt("N=1000000_fracsize_vs_t.txt"))
labels.append(r'$N=100000$ Achlioptas, $m=2$')
colors.append(tableau20[6])

'''
plot_data.append(genfromtxt("N=1000_fracsize_vs_t.txt"))
labels.append(r'$N=1000$')
colors.append(tableau20[0])

plot_data.append(genfromtxt("N=10000_fracsize_vs_t.txt"))
labels.append(r'$N=10000$ (Vieja)')
colors.append(tableau20[2])

plot_data.append(genfromtxt("N=100000_fracsize_vs_t.txt"))
labels.append(r'$N=100000$')
colors.append(tableau20[4])

plot_data.append(genfromtxt("N=100_fracsize_vs_t.txt"))
labels.append(r'$N=100$')
colors.append(tableau20[6])
'''

fig = plt.figure()
i=0
plt.plot((plot_data[i][:,0]), plot_data[i][:,1], '--',#, plot_data[i][:,2],
                label = labels[i],color = 'k', linewidth=2.0,
            ) 
i=1
plt.plot((plot_data[i][:,0]), plot_data[i][:,1], #marker='o', facecolors='none',# plot_data[i][:,2],
                label = labels[i],color = colors[i], linewidth=2.0#, s=40.0,
            ) 

'''
for i in range(0,len(plot_data)):
    plt.plot((plot_data[i][:,0]), plot_data[i][:,1], '--',#, plot_data[i][:,2],
                label = labels[i],color = colors[i], linewidth=2.0,
            )    
'''
plt.legend()

#plt.yscale('log')
# Formating:
plt.xlim(0.0,1.0)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\vert C \vert/N$')
plt.ylim(0,1.1)

plt.gcf().subplots_adjust(bottom=0.15)
plt.show()
#plt.grid(True)
#plt.rc('grid', linestyle='--', color='gray')
plt.savefig("perc.pdf")