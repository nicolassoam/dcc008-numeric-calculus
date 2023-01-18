import numpy as np
import matplotlib.pyplot as plt
import math as mh
from timeit import default_timer as timer


# Define tipo da variavel
TYPE = np.float128

# Método de Newton para a solução de Sistemas Não Lineares

# Recebe como entrada a matriz de sistema de funções não lineares
# A matriz jacobiana associada a essas funções
# A aproximação inicial
# Número máximo de iterações

def newtonRaphson(a,b,x,e,N,max,alpha):
  # Funcao F(xi)
  func = np.zeros(2*N,dtype=TYPE)

  # Matriz Jacobiana
  J = np.zeros((2*N,2*N),dtype=TYPE)

  # Integrais aproximadas por Newton Cotes
  g = newtonCotes(f,a,b,N,1500)

  # Vetor auxiliar
  s = np.zeros(2*N, dtype=TYPE)

  i = 0;

  while(np.linalg.norm(x,np.inf)>e):
    func = calculaF(x[0:N],x[N:2*N],N,g)
    J = jacobiana(x,e,N,a,b,g,alpha)
    s = decLU(J,-func,2*N)

    # Posicoes 0 a N do vetor x possuem os pesos w
    x[0:N] = s[0:N] + x[0:N]

    # Posicoes N a 2N do vetor x possuem os pesos t
    x[N:2*N] = s[N:2*N] + x[N:2*N]

    i+=1
    if(i == max):
      print("Atingiu o número máximo de iterações.")
      return x

  return x


# Define um chute inicial
# Dependendo dos valores, não funciona corretamente
def chuteInicial(a,b,N):
  x0 = np.zeros(2*N,dtype=TYPE)
  w0 = np.zeros(N, dtype=TYPE)
  t0 = np.zeros(N, dtype=TYPE)

  # Calcula w0
  for i in range(int(N/2)):
      w0[i] = ((b-a)/(2*N))*(i+1)
      x0[i] = w0[i]
      t0[i] = a+ (i+1)*(w0[i]/2)
      x0[i+N] = t0[i]
      
      w0[N-1-i] = w0[i]
      x0[N-1-i] = w0[N-1-i]
      t0[N - 1 - i] = (a+b) - t0[i]
      x0[N-1-i + N] = t0[N-1-i]


  if(N%2 == 0):
    return x0,w0,t0
  else:
    t0[int(N/2)] = (a+b)/2
    x0[int(N/2)+N] = t0[int(N/2)]
    w0[int(N/2)] = (a-b)/2
    x0[int(N/2)] = w0[int(N/2)]
    return x0,w0,t0

# Calcula integrais polinomias sem utilizar newton cotes
def calculaIntegrais(a,b,N):
  s = np.zeros(2*N,dtype=TYPE)
  s[0] = b-a
  for i in range(1,2*N):
    s[i] = ((b**(i+1))/(i+1)) - ((a**(i+1))/(i+1))
  return s

# Newton Cotes (Simpsons 3/8)
def newtonCotes(function, lower_bound, upper_bound,N, intervals):
  g = np.zeros(2*N,dtype=TYPE)

  for j in range(0,2*N):
    I = 0.
    # calculating step size
    h = (b - a) / intervals
    
    # Finding sum 
    I = f(a,j) + f(b,j)
    
    for i in range(1,intervals):
        k = a + (i*h)
        
        if i%3 == 0:
            I = I + (2 * f(k,j))
        else:
            I = I + (3 * f(k,j))
    
    # Finding final integration value
    I = I * (3 * (h / 8))
    g[j] = I


  # Return the result of the integration
  return g
  
# expoente para newton cotes
def f(n,i):
  return n**i

# Calcula o vetor F
def calculaF(w0,t0,N,integral):
  func = np.zeros(2*N,dtype=TYPE)
  sum = 0.
  for j in range(2*N):
    sum = 0.
    for i in range(N):
      sum  += (w0[i]*(t0[i]**j))

    func[j] = sum - integral[j]
  
  return func

# Calcula a Matriz Jacobiana
def jacobiana(x,eps,N,a,b,integral,alpha):
  jac = np.zeros((2*N,2*N),dtype=TYPE)

  for i in range(2*N):
    x_ = x.copy()
    x_[i] = x_[i] + eps;
    for j in range(2*N):
      f1_ = calculaF(x_[0:N],x_[N:2*N],N,integral)
      f1 = calculaF(x[0:N],x[N:2*N],N,integral)
      jac[j,i] = (f1_[j]-f1[j])/eps
    x_[i] = x[i]
  for i in range(2*N):
    jac[i,i] +=alpha

  return jac


#Gauss com pivoteamento parcial
def pivoteamentoParcial(A,b,n):
  aux = np.zeros(n,dtype=TYPE);
  r = -1;
  m = 0

  for k in range(0,n-1):
    w = np.abs(A[k][k]);
    for j in range(k,n):
      if(np.abs(A[j][k])>w):
        w = np.abs(A[j][k]);
        r = j

    if (r!=-1): 
      A[[r,k]] = A[[k,r]];
      x = b[r]
      b[r] = b[k]
      b[k] = x

    m = A[k+1:,k]/A[k,k];
    A[k+1:,k:n] = A[k+1:,k:n] - A[k,k:n]*m[:,np.newaxis];
    b[k+1:] = b[k+1:] - m*b[k];

  aux = retroSub(A,b,n);

  return aux;

# Gauss sem pivoteamento
def eliGauss(A,b,n):
  
    aux = np.zeros(n,dtype=TYPE);
    m = 0

    for k in range(0,n):
        m = A[k+1:,k]/A[k,k];
        A[k+1:] = A[k+1:] - A[k]*m[:,np.newaxis];
        b[k+1:] = b[k+1:] - m*b[k];
    aux = retroSub(A,b,n);
   
    return aux

#Decomposição LU
def decLU(A,b,n):
  m = 0
  #Matrizes Extra
  L = np.eye(n,dtype=TYPE);
  U = np.zeros((n,n),dtype=TYPE);


  for k in range(0,n):
    m = A[k+1:,k]/A[k,k];
    L[k+1:,k] = m
    A[k+1:] = A[k+1:] - A[k]*m[:,np.newaxis];
    
  U = A;
  sub = subst(L,b,n)
  retrosub = retroSub(U,sub,n)

  return retrosub

# Substituição
def subst(A,b,n):
  aux = np.zeros(n,dtype=np.float64());
  aux[0] = b[0]/A[0][0];
  for i in range(1,n):
    s = b[i];
    for j in range(0,i):
      s = s - A[i][j]*aux[j];
    aux[i] = s/A[i][i];
  return aux;

# Retro Substituição
def retroSub(A,b,n):
    aux = np.zeros(n,dtype=TYPE);
    aux[n-1] = b[n-1]/A[n-1][n-1];

    for i in range(n-2,-1,-1):
      s = b[i];
      for j in range(i+1,n):
        s = s - A[i][j]*aux[j];
      aux[i] = s/A[i][i];
    return aux;

# Função exata
def exata(a,b,x,expOr):
  if(expOr):
    return np.exp(a*x+b)
  else:
    return np.sqrt((a**2+b))

# Gauss Legendre
def GaussLegendre(a,b,N,eps,alpha,max):
  start = timer()
  x0 = chuteInicial(a,b,N)
  weights = newtonRaphson(a,b,x0[0].copy(),eps,N,max,alpha)
  sum = 0.
  print("O chute inicial foi : ")
  print(x0[0])
  print("")
  print("Os pesos obtidos foram ([0,N] = w, [N, 2N] = t): ")
  print(weights)
  print("")

  for i in range(N):
    w = weights[i]
    t = weights[i+N]
    sum += w* exata(a,b,t,False) 
  end = timer()
  tempo = end - start
  print("Tempo de execução: " + str(tempo))
  return sum

def integralExata(a,b):
  e = (np.exp((a*b)+b) - np.exp((a**2)+b))/a
  return e

def integralExata2(a,b):
  root = (np.sqrt((a**2+b))*b) - (np.sqrt((a**2+b))*a)
  return root

def erro(exata,x):
  return np.abs(np.subtract(exata,x))/x
