import scipy.optimize as opt

def promedio_dif (y):
    k_avg = 6
    NODE_NR = 300
    K_MIN = 2
    K_MAX = NODE_NR
    
    suma = 0
    for k in range(K_MIN,K_MAX):
        suma = suma + k**(-y)
        
    suma2 = 0
    for k in range(K_MIN,K_MAX):
        suma2 = suma2 + k* k**(-y)
        
    return k_avg - suma2 / suma, suma, suma2

def promedio_root (y):
    return promedio_dif(y)[0]



ini_gamma = 2.5
ROOT = opt.root(promedio_root, ini_gamma, method='hybr')

print(ROOT.x)