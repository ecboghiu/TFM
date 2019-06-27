import matplotlib.pyplot as plt
import networkx as nx
import collections
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

#sns.set()
#sns.set_context('talk')
#sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 1.5})


import numpy as np
from pylab import genfromtxt


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


def chop_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#This cmap is applied in to node coloring 
cmap_aux = LinearSegmentedColormap.from_list("", ["red","yellow"])

#In gray scale, cyan becomes almost white, this is why we chop the begining of the color map
cyan_purple_cmap = chop_colormap(cmap_aux, 0.05, 1)

def spectrum_to_color(w):
    colors_cyan_purple = cyan_purple_cmap(np.linspace(0, 1, 1001))
    w_min = min(w)
    w_max = max(w)
    return [colors_cyan_purple[int(1000*(w[i] - w_min)/(w_max - w_min))] for i in range(len(w))]



'''
import seaborn as sns
sns.set()
#sns.set_context('talk')
sns.set_context("notebook", font_scale=1.1, rc={"lines.linewidth": 1.5})
'''


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

plot_data, labels, colors, stat_hist = [], [], [], []
data_labels_colors = [plot_data, labels, colors, stat_hist]

def hist_3d (filename_in, filename_out):
    datos = (np.genfromtxt(filename_in,
                           skip_header=5))
    xdat=datos[:,1]
    ydat=datos[:,2]
    #NODE_NR = 500
    G=nx.Graph()
    #for i in range(0,NODE_NR):
    #    G.add_node(i)
    '''
    options = {
        'node_color': 'black',
        'node_size': 50,
        'line_color': 'grey',
        'linewidths': 0,
        'width': 0.1,
    }
    '''
    hist_surf = []
    for x, y in zip(xdat,ydat):
        G.add_edge(x,y)
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = zip(*degreeCount.items())
        
        np_deg = np.array(deg)
        np_cnt = np.array(cnt)
        norm = np.sum(np_cnt)   
        np_prob = np_cnt/norm
        average = np.sum(np_deg*np_prob)
        k2 = np.sum(np_deg*np_deg*np_prob)
        
        hist_surf.append([deg, cnt, average, k2])
        
    #fig = plt.figure()
    #i=10
    #plt.xlabel(r"$\mathrm{log} k$")
    #plt.ylabel(r"$\mathrm{log} f$")
    #for j in range(0,len(hist_surf)): 
    #j=(len(hist_surf)-1)
    #plt.plot((hist_surf[j][0][:]),(hist_surf[j][1][:]))
    

    #aux = np.linspace(0,len(hist_surf)/NODE_NR,len(hist_surf),dtype=np.float32)
    #aux2 = np.ones(len(hist_surf))
    
    #ax = fig.add_subplot(111, projection='3d')
    #ax = Axes3D(fig)
    '''
    file = open(filename_out,"w")
    for i in range(0,len(hist_surf)):
        #if (i % (len(hist_surf)/3) == 0 ):
        for j in range(0,len(hist_surf[i][0])):   
            if (j % (len(hist_surf[i][0])/(len(hist_surf[i][0])/10)) == 0 ):
                file.write(str(aux[i])+" "+str(np.log(hist_surf[i][0][j]))+" "+str(np.log(hist_surf[i][1][j]))+"\n")
    file.close()  
    '''    
    
    #print("End of function"+filename_in+filename_out)
    return hist_surf

def take_from_txt (out, filename, filename_hist, color ):
    out[0].append(genfromtxt(filename,  skip_header=5))
    out[1].append(['r_'+filename,'p_'+filename])
    color_code = color
    out[2].append([tableau20[color_code], tableau20[color_code+1],tableau20[color_code]])
    
    
    file_in = filename_hist
    file_out = "../data/hist3D_"+filename
    out[3].append(hist_3d(file_in,file_out))
    
    
   
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


#variation with h
#take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_a=5_sig=0.08_h=.1.txt", 6)
#take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_a=5_sig=0.08_h=.05.txt", 4)
#take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_a=5_sig=0.08_h=1..txt", 0)
#take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_a=5_sig=0.08_h=0.01.txt", 2)
lambas=['0.1','0.09','0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']
lambas=['0.04','0.03','0.02']
'''
lambas=['0.07','0.06','0.058']
lmmm = [str(0.01+0.002*i).format("%.3d") for i in range(0,45)]

lmmm
lambas= ['0.01',
 '0.012',
 '0.014',
 '0.016',
 '0.018',
 '0.02',
 '0.022',
 '0.024',
 '0.026',
 '0.028',
 '0.03',
 '0.032',
 '0.034',
 '0.036',
 '0.038',
 '0.04',
 '0.042',
 '0.044',
 '0.046',
 '0.048',
 '0.05',
 '0.052',
 '0.054',
 '0.056',
 '0.058',
 '0.06',
 '0.062',
 '0.064',
 '0.066',
 '0.068',
 '0.07',
 '0.072',
 '0.074',
 '0.076',
 '0.078',
 '0.08',
 '0.082',
 '0.084',
 '0.086',
 '0.088',
 '0.09',
 '0.092',
 '0.094',
 '0.096',
 '0.098',
 '0.1']
'''

#lambas_list=['1','0.5','0.1','0.08','0.05']
lam = lambas
lambas=lambas
#alpha_list=['-20','-10','-1','0','1','5','20']
alpha_list=['-20','-5','0','5','20']
def data_section_lam(lam, alpha, data_labels_colors):
    for i in range(len(alpha)):
        take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=2000_m=1_a="+(alpha[i])+"_sig="+lam+".txt",
                  "../data.19.6.11.achlioptas/FG_N=2000_m=1_a="+(alpha[i])+"_sig="+lam+"_EDGELIST.txt", 2*i)
def data_section_alph(lam, alpha, data_labels_colors):
    for i in range(len(lam)):        
        take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a="+(alpha)+"_sig="+(lam[i])+".txt",
                  "../data/FG_N=2000_m=1_a="+(alpha)+"_sig="+(lam[i])+"_EDGELIST.txt", 2*i)

#data_section_lam('1', alpha, data_labels_colors)
        
alpha_number = 3
data_section_alph(lam, alpha_list[alpha_number], data_labels_colors)
#data_section_lam('0.5',alpha_list,data_labels_colors)
        
#estas son las chulas
'''
..

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.5.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.5_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.1.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.1_EDGELIST.txt", 0)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.09.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.09_EDGELIST.txt", 2)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.08.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.08_EDGELIST.txt", 4)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.07.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.07_EDGELIST.txt", 6)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.06.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.06_EDGELIST.txt", 8)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.05.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.05_EDGELIST.txt", 10)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.04.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.04_EDGELIST.txt", 12)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.03.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.03_EDGELIST.txt", 14)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.02.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.02_EDGELIST.txt", 16)

take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.01.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.01_EDGELIST.txt", 18)
'''
'''
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.009.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.009_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.008.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.008_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.007.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.007_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.006.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.006_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.005.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.005_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.004.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.004_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.003.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.003_EDGELIST.txt", 0)

...
'''

# estos son las de 2000 pero sin histeresis porque se quedaron bloqueadas
'''
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.01.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.01_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.009.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.009_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.008.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.008_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.007.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.007_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.006.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.006_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.005.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.005_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.004.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.004_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.003.txt",
                  "../data/FG_N=2000_m=1_a=5_sig=0.003_EDGELIST.txt", 0)
'''

# estos son para ell super alto
'''
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.05.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.05_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.01.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.01_EDGELIST.txt", 0)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.007.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.007_EDGELIST.txt", 2)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.004.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.004_EDGELIST.txt", 4)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.003.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.003_EDGELIST.txt", 6)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.0025.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.0025_EDGELIST.txt", 8)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.002.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.002_EDGELIST.txt", 10)
take_from_txt(data_labels_colors, "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.001.txt",
                  "../data.19.6.11.achlioptas/FG_N=500_m=1_a=5_sig=0.001_EDGELIST.txt", 12)
'''


'''
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.07.txt", 2)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.06.txt", 4)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.05.txt", 6)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.04.txt", 8)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.03.txt", 10)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.02.txt", 12)
take_from_txt(data_labels_colors, "../data/FG_N=2000_m=1_4_a=5_sig=0.01.txt", 14)
'''
#take_from_txt(data_labels_colors, "../data/ach_change_with_k/FG_N=1000_m=1_2_a=0_sig=0.08.txt", 6)

#filename = "../EPES_N=1000_tribe_sig=5.txt"
#data_labels_colors[0].append(genfromtxt(filename,  skip_header=5))
#data_labels_colors[1].append(['r_'+filename,'p_'+filename])
#color_code = 0
#data_labels_colors[2].append([tableau20[color_code],
#                    tableau20[color_code+1],tableau20[color_code]]) 


# DOWNSAMPLING SYN
'''
NUMBER_OF_SAMPLES = 1000
loop_nr=len(plot_data)
for i in range(0,loop_nr):
    #plot_data[i] = plot_data[i][::(int(len(plot_data[i])/NUMBER_OF_SAMPLES))].copy()
    plot_data[i] = plot_data[i][::100].copy()
'''

def r_to_color(r_list):
    #colors_cyan_purple = cyan_purple_cmap(np.linspace(0, 1, 1001))
    aux_len = len(r_list)
    colors_aux = cm.rainbow(np.linspace(0, 1, aux_len))
    return [colors_aux[int(aux_len*r_list[i])] for i in range(0,aux_len)]


def plot_epes(data_labels_colors): 
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9,5))
    
    
    histograms = data_labels_colors[3]
    
    NODE_NR = 2000
    
    
    '''
    j_0=int(2.5*NODE_NR)
    j_1=int(3.5*NODE_NR)
    for j in [j_0,j_1]:#range(j_0, j_1):
        ax[1][0].plot((histograms[0][j][0][:]),(histograms[0][j][1][:]),
                      label = "$\ell$="+str(float(j/NODE_NR)))
    ax[1][0].legend()
	'''
    
    
    ax[1].set_xlabel(r'$\ell$')
    ax[1].set_ylabel(r'$\alpha$')
    #color_aux = np.linspace(0, 1, len(plot_data[0][:,0]))
    #coloress = spectrum_to_color(plot_data[0][:,2])
    
    #lambas=['0.5','0.1','0.09','0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']
    '''
    lambas=['0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']
    for i in range(0,len(plot_data)):
        temp_len = int(len(plot_data[i][:,0])/2)
        #coloress = r_to_color(plot_data[i][:temp_len,2])
        #for j in range(0,int(len(plot_data[i][:,0]/2))):
        aux_cb = ax[1][1].scatter(plot_data[i][:temp_len,0],
                    float(lambas[i])*np.ones(temp_len),
                    #linewidth=0.5, #marker='+',
                    s = 20,
                    label = labels[i][0],
                    c = plot_data[i][:temp_len,2],
                    cmap=plt.cm.get_cmap('viridis'))
    
    #ax[1][1].set_xscale('log')
    
    ax[1][1].set_xlim(0.1 , 30.0)
    ax[1][1].set_ylim(0.0 ,  1.0)
    cb_aux= plt.colorbar(aux_cb, ax=ax[1][1])
    cb_aux.set_label(r'$r$')
    #sm = plt.cm.ScalarMappable(cmap=cm.rainbow(np.linspace(0, 1, 1000)), norm=plt.Normalize(vmin=0, vmax=1))
    #sm._A = []
    #cb1 = plt.colorbar(aux)
    #cb1.set_label(r'$r$')
    '''
    ax[0].set_xlabel(r'$\ell$')
    ax[0].set_ylabel(r'$r,\kappa$')
    
    for i in range(0,len(plot_data)):
        temp_len = int(len(plot_data[i][:,0])/2)
        ax[0].plot((plot_data[i][:temp_len,0]), plot_data[i][:temp_len,2],# np.array(plot_data[i][:,3])/np.sqrt(250),
                    linewidth=1,color = colors[i][0], #marker='+',
                    label=r'$r$ ($\sigma=$'+lambas[i]+")")#labels[i][0])#,color = colors[i][0])
        ax[1].plot((plot_data[i][:temp_len,0]), plot_data[i][:temp_len,1],
          color = colors[i][0],
          label = r'$\sigma=$'+lambas[i])#labels[i][1])#, color = colors[i][1], linewidth=1)
        
        aux_array = np.array(histograms[i])
        aux_temp_len = len(aux_array)
        k2byk = (aux_array[:,3])/(aux_array[:,2])
        ax[1].plot((plot_data[i][:aux_temp_len,0]), k2byk, color = colors[i][0],
          label=r'$r$ ($\sigma=$'+lambas[i]+")")
        ax[0].plot((plot_data[i][:aux_temp_len,0]), k2byk/50, '--', color = colors[i][0],#(1/k2byk)*2/np.pi/1.0/75,
          label=r'$\frac{\kappa}{50}$ ($\sigma=$'+lambas[i]+")") 
        
    
    #ax[0].set_xscale('log')
    
    #ax[0].set_xlim(1,30)
    
    ax[0].legend().set_draggable('True')
    ax[1].legend().set_draggable('True')
    
    #plt.yscale('log')
    # Formating:
    #plt.xlim(0,1.1)
    '''
    ax[0][0].set_xlabel(r'($\ell$)')
    ax[0][1].set_xlabel(r'($\ell$)')
    ax[0][1].set_ylabel(r'$\vert C \vert /N$')
    ax[0][0].set_ylabel(r'r')
    '''
    #ax[0][0].set_ylim(0,1.1)
    #ax[0][1].set_ylim(0,1.1)
    #plt.title(r'Synchronization and maximum component fraction')
    
    plt.tight_layout()
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig("r_dep.pdf")
    plt.show()
    
def plot_epes_222(data_labels_colors): 
    fig, ax = plt.subplots(nrows=2, ncols=len(alpha_list), figsize=(20,8))
    
    NODE_NR = 2000

    col=0
    for alfa_aux in alpha_list:
        data_labels_colors.clear()
        plot_data, labels, colors, stat_hist = [], [], [], []
        data_labels_colors = [plot_data, labels, colors, stat_hist]
        data_section_alph(lambas_list, alfa_aux, data_labels_colors)
        ax[0][col].set_xlabel(r'$\ell$')
        ax[0][col].set_ylabel(r'$r$')
        ax[1][col].set_xlabel(r'$\ell$')
        ax[1][col].set_ylabel(r'$\vert C \vert /N$')
        
        for i in range(0,len(plot_data)):
            ax[0][col].plot((plot_data[i][:,0]), plot_data[i][:,2],# np.array(plot_data[i][:,3])/np.sqrt(250),
                    linewidth=1, #marker='+',
                    label = labels[i][0])#,color = colors[i][0])
            ax[1][col].plot((plot_data[i][:,0]), plot_data[i][:,1],
          label = labels[i][1])#, color = colors[i][1], linewidth=1)
        
        
        col = col + 1
    
    
    ax[0][0].legend().set_draggable('True')
    ax[0][1].legend().set_draggable('True')
    
    #
    plt.tight_layout()
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig("r_dep_222.pdf")
    plt.show()
    
def plot_weff_corr():
    datos = (genfromtxt("../data/FG_N=2000_m=1_a=10_sig=0.08_EDGELIST.txt",
                       skip_header=1))
    
    plt.xlabel(r"$\omega_{{eff}}$");
    plt.ylabel(r"$\omega_{{eff}}$");
    xdat=datos[:,2]
    ydat=datos[:,3]
    
    #colors = cm.rainbow(np.linspace(0, 1, len(ydat)))
    #for x, y, c in zip(xdat,ydat,colors):
    #    plt.scatter(x, y, color=c, marker='+')
    sns.jointplot(x=x, y=y, kind='hex')
    plt.show()
    
def plot_r_thing(data_labels_colors):
    plt.figure(figsize=(2*one_column_figure_size * golden_ration/1.4, 2*one_column_figure_size))
    #lambas1=['0.10','0.09','0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']
    #lambas=['0.05','0.01','0.007','0.004','0.003','0.0025','0.002','0.001']
    #lambas=['1.00','0.50','0.10','.09','0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']
    lambas=lambas_list
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$r$")
    plt.title(r"$\alpha=5$")
    plt.xlim(0,4)
    plt.ylim(0,1.1)
    #plt.xscale('log')
    for i in range(0,len(plot_data)):
        plt.plot((plot_data[i][:,0]), plot_data[i][:,2],# np.array(plot_data[i][:,3])/np.sqrt(250),
                linewidth=1, #marker='+', markeredgewidth=0.5,
                label = r"$\sigma=$"+str(lambas[i]))#labels[i][0], #color = colors[i][0])
    plt.tight_layout()
    plt.legend(loc='upper left').set_draggable('True')
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig("graph_aux_r.pdf")
    
def plot_perc_thing(data_labels_colors):
    plt.figure(figsize=(2*one_column_figure_size * golden_ration/1.4, 2*one_column_figure_size))
    #lambas=['0.10','0.09','0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']
    lambas=lambas_list
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\vert C \vert /N$")
    plt.title(r"Con proceso de Achlioptas")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlim(0,1)
    plt.ylim(0,1.1)
    aux_len=int(len(plot_data[0][:,0])/2)
    for i in range(0,len(plot_data)):
        plt.plot((plot_data[i][:aux_len,0]), plot_data[i][:aux_len,1],
                 linewidth=1,
                label = r"$\alpha=$"+str(alpha_list[i]))#label = r"$\sigma=$"+str(lambas[i]))#, color = colors[i][1], linewidth=1)
    plt.tight_layout()
    plt.legend(loc='lower right').set_draggable('True')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig("graph_aux_perc.pdf")
    
'''
def plot_hist_thing(data_labels_colors):
    histograms = data_labels_colors[3]
    plt.xlabel(r"$k$")
    plt.ylabel(r"$\propto P(k)$")
    #plt.yscale('log')
    #plt.xscale('log')
    NODE_NR = 2000
    j_0=int(8.8*NODE_NR)
    j_1=int(9.0*NODE_NR)
    j_2=int(11.6*NODE_NR)
    j_3=int(11.7*NODE_NR)
    rango = [j_0,j_1,j_2,j_3]
    for j in rango:#range(j_0, j_1):
        plt.plot((histograms[0][j][0][:]),(histograms[0][j][1][:]),
                      label = "$\ell$="+str(float(j/NODE_NR)))
    plt.legend()
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig("graph_aux_hist.pdf")
'''

def plot_phase_diag(data_labels_colors):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\sigma$')
    #color_aux = np.linspace(0, 1, len(plot_data[0][:,0]))
    #coloress = spectrum_to_color(plot_data[0][:,2])

    NODE_NR = 2000
    
    #lambas=['0.5','0.1','0.09','0.08','0.07','0.06','0.05','0.04','0.03','0.02','0.01']#,'0.009','0.008','0.007','0.006','0.005','0.004','0.003']
    #t1t2 = [[0,0], [2.205,2.2608], [2.573,2.657],
    #        [2.991,3.233], [3.593,4.052], [4.497,5.294], [6.00,7.420],
    #        [8.97,11.642], [17.679,23.689]]
    '''
    t1t2 = []
    for i in range(0,len(plot_data)):
        aux_a = plot_data[i][:,2] # sinc
        aux_b = plot_data[i][:,0] # ell
        #t2  = float(aux_b[np.argwhere(aux_a>0.3)[0]])
        #t1  = float(aux_b[np.argwhere(aux_a>0.3)[-1]])
        t2  = int(np.argwhere(aux_a>0.3)[0])
        t1  = int(np.argwhere(aux_a>0.3)[-1])
        #if ((aux_b[t2]-aux_b[t1])*NODE_NR)>100:
        if t1 < t2:
            t1t2.append([t1,t2])
        if t1 > t2:
            t1t2.append([t2,t1])
    '''
        
    for i in range(0,len(plot_data)):
        temp_len = int(len(plot_data[i][:,0])/2)
        #coloress = r_to_color(plot_data[i][:temp_len,2])
        #for j in range(0,int(len(plot_data[i][:,0]/2))):
        #aux_cb = ax.scatter(plot_data[i][:temp_len,0],
        #            float(lambas[i])*np.ones(temp_len), marker='s', #linewidth=0.5, #marker='+',
        #            s = 20, #label = labels[i][0],
        #            c = plot_data[i][:temp_len,2],
        #            cmap=plt.cm.get_cmap('plasma'),
        #            norm=plt.Normalize(vmin=0, vmax=1))
        
        plt_len = len(plot_data[i][temp_len:,0])
        aux_cb=ax.scatter(plot_data[i][temp_len:,0],
            float(lambas[i])*np.ones(plt_len), marker='s', #linewidth=0.5, #marker='+',
            s = 20, #label = labels[i][0],
            c = plot_data[i][temp_len:,2],
            cmap=plt.cm.get_cmap('viridis'),
            norm=plt.Normalize(vmin=0, vmax=1))

        
        '''
        #t1 = int(t1t2[i][0]*NODE_NR)
        #t2 = int(t1t2[i][1]*NODE_NR)
        t1 = t1t2[i][0]
        t2 = t1t2[i][1]
        ax.scatter(plot_data[i][t1:t2,0],
                    float(lambas[i])*np.ones(abs(t2-t1)), marker='s',
                    #linewidth=0.5, #marker='+',
                    s = 10, #label = labels[i][0],
                    c = 'cyan')
        '''
        
    #sm = plt.cm.ScalarMappable(cmap=cm.rainbow(np.linspace(0, 1, 1000)), norm=plt.Normalize(vmin=0, vmax=1))
    #sm._A = []
    #cb1 = plt.colorbar(aux)
    #cb1.set_label(r'$r$')

    
    #ax[1][1].set_xscale('log')
    #plt.yscale('log')
    #plt.xscale('log')
    ax.set_xlim(0.1 , 25.0)
    ax.set_ylim(0.01 ,  0.1)
    cb_aux= plt.colorbar(aux_cb, ax=ax)
    cb_aux.set_label(r'$r$')
    #plt.savefig("phase_diag.pdf")
    plt.tight_layout()
    plt.savefig('destination_path.png', format='png', dpi=1200)
    
def plot_both_thing(data_labels_colors):
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(5,4))
    
    common_x_lim = 2
    
    lambas=lambas_list
    #ax[0].set_xlabel(r"$\ell$")
    ax[0].set_ylabel(r"$r$")
    ax[0].set_title(r"$\alpha="+alpha_list[alpha_number]+"$")
    #ax[0].set_xlim(0,common_x_lim)
    ax[0].set_ylim(0,1)
    ax[0].axes.xaxis.set_ticklabels([])
    #plt.xscale('log')
    for i in range(0,len(plot_data)):
        ax[0].plot((plot_data[i][:,0]), plot_data[i][:,2],# np.array(plot_data[i][:,3])/np.sqrt(250),
                linewidth=1, #marker='+', markeredgewidth=0.5,
                label = r"$\sigma=$"+str(lambas[i]))#labels[i][0], #color = colors[i][0])
    #ax[0].tight_layout()
    #ax[0].legend().set_draggable('True')
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.savefig("graph_aux_r.pdf")

    
    ax[1].set_xlabel(r"$\ell$")
    ax[1].set_ylabel(r"$\vert C \vert /N$")
    #plt.yscale('log')
    #plt.xscale('log')
    ax[1].set_xlim(0,common_x_lim)
    ax[1].set_ylim(0,1)
    
    for i in range(0,len(plot_data)):
        aux_len=int(len(plot_data[i][:,0])/2)
        ax[1].plot((plot_data[i][:aux_len,0]), plot_data[i][:aux_len,1],
                 linewidth=1,
                label = r"$\sigma=$"+str(lambas[i]))#, color = colors[i][1], linewidth=1)
    #ax[1].tight_layout()
    ax[1].legend().set_draggable('True')
    #ax[1].gcf().subplots_adjust(bottom=0.15)
    #ax[1].savefig("graph_aux_perc.pdf")
    
    plt.tight_layout()
    
    
    
#plot_both_thing(data_labels_colors)

#plot_r_thing(data_labels_colors)
#plot_perc_thing(data_labels_colors)

#plot_phase_diag(data_labels_colors)
plot_epes(data_labels_colors)
#p
#plot_weff_corr()

