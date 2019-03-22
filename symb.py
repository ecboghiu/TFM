import scipy.optimize as opt
import scipy as sc
from scipy.stats import norm

def promedio_SF (y):
    k_avg = 6
    NODE_NR = 1000
    K_MIN = 2
    K_MAX = NODE_NR
    
    suma = 0
    for k in range(K_MIN,K_MAX):
        suma = suma + k**(-y)
        
    suma2 = 0
    for k in range(K_MIN,K_MAX):
        suma2 = suma2 + k * k**(-y)
        
    suma3 = 0
    for k in range(K_MIN,K_MAX):
        suma3 = suma3 + k * k * k**(-y)
        
    return k_avg - suma2 / suma, suma, suma2/suma, suma3/suma

def promedio_root (y):
    return promedio_SF(y)[0]

ini_gamma = 2.5
ROOT = opt.root(promedio_root, ini_gamma, method='hybr')

K_meanfield = 2.0/sc.pi/norm.pdf(0)

avg_k  = promedio_SF(ROOT.x)[2]
avg_k2 = promedio_SF(ROOT.x)[3]
print("K_SF=", K_meanfield * avg_k / avg_k2/1000)

print("K_ER=", K_meanfield/1000)


print("gamma=",ROOT.x)