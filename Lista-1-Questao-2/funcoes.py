import numpy as np
import matplotlib.pyplot as plt

#Solução exata

def calc_c2(eps, h):
  return (np.exp((-1/np.sqrt(eps)))-1)/(np.exp((1/np.sqrt(eps))) - (np.exp((-1/np.sqrt(eps)))-1));

def calc_c1(c2):
  return (-1 - c2);

def exact_solution(x,eps,h):
  c2 = calc_c2(eps, h);
  c1 = calc_c1(c2);
  e_u = c1*np.exp((-x/np.sqrt(eps))) + c2*np.exp((x/np.sqrt(eps)))+1;
  return e_u;

t = np.linspace(0,1,50);
h = 0.020;

def TDMA(a,b,c,d,ordem,eps,h):
    
    a[:ordem] = -eps;
    c[:ordem] = -eps;
    
    b[:ordem+1] = (2*eps+(h**2));
    d[:ordem+1] = (h**2);

    # for n in range(0,ordem):
    #   a[n] = -eps;
    #   c[n] = -eps;
    # for n in range(0,ordem+1):
    #   b[n] = (2*eps+(h**2));
    #   d[n] = (h**2);
    # d[0] = 0;
    # d[ordem] = 0;

    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
      
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    xc[0] = 0.;
    xc[ordem] = 0.;

    return xc

#diagonal inferior
a = np.zeros(49);

#diagonal central
b = np.zeros(50);

#diagonal superior
c = np.zeros(49);

#vetor solucao
d = np.zeros(50);