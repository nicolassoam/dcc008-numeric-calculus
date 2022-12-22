import numpy as np;
import matplotlib.pyplot as plt
from timeit import default_timer as timer

class SistLinear:

  def __init__(self,n=None):
    if(n is None):
      n = 10
    self.n = n;
    #Matriz nxn
    self.A = np.zeros((n,n), dtype=np.float64())
    #Vetor solucao
    self.b = np.zeros(n,dtype=np.float64());

    #Gera matriz
    for i in range(n):
      self.b[i] = 1/(i+n+1);
      for j in range(n):
        self.A[i][j] = 1/(j+i+1);

    #Determinante de A
    detA = 0
    detA = np.linalg.det(self.A);
    print("det = " + str(detA));

    #Verifica se a solução é única
    if (detA == 0): return;
  
  def retornaMatriz(self):
    return self.A;
  
  def retornaSolucao(self):
    return self.b;

  #Retro-Substituição
  def retroSub(self,A,b,n):
    #A,b,n = self.A.copy(),self.b.copy(),self.n
    aux = np.zeros(n,dtype=np.float64());
    aux[n-1] = b[n-1]/A[n-1][n-1];

    for i in range(n-2,-1,-1):
      s = b[i];
      for j in range(i+1,n):
        s = s - A[i][j]*aux[j];
      aux[i] = s/A[i][i];
    return aux;
  
  #Substituição
  def subst(self,A,b,n):
    #A,b,n = self.A.copy(),self.b.copy(),self.n
    aux = np.zeros(n,dtype=np.float64());
    aux[0] = b[0]/A[0][0];
    for i in range(1,n):
      s = b[i];
      for j in range(0,i):
        s = s - A[i][j]*aux[j];
      aux[i] = s/A[i][i];
    return aux;

  #Eliminação de Gauss
  def eliGauss(self,A=None,b=None,n=None):
    A,b,n = self.__verifyParameters(A,b,n)
    aux = np.zeros(n,dtype=np.float64());
    m = 0

    start = timer()

    # for k in range(0,n):
    #   for i in range(k+1,n):
    #     m = A[i][k]/A[k][k];
    #     for j in range(k,n):
    #       A[i][j] = A[i][j] - A[k][j]*m;
    #     b[i] = b[i] - m*b[k];
   
    for k in range(0,n):
        m = A[k+1:,k]/A[k,k];
        A[k+1:] = A[k+1:] - A[k]*m[:,np.newaxis];
        b[k+1:] = b[k+1:] - m*b[k];
    aux = self.retroSub(A,b,n);

    end = timer()

    return aux, self.erro(A,b,aux),end-start;

  #Gauss com pivoteamento parcial
  def pivoteamentoParcial(self,A=None,b=None,n=None):
    A,b,n = self.__verifyParameters(A,b,n)
    aux = np.zeros(n,dtype=np.float64());
    r = -1;
    m = 0

    start = timer()

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

      # for i in range(k+1,n):
      #     m = A[i][k]/A[k][k];
      #     for j in range(k,n):
      #       A[i][j] = A[i][j] - m*A[k][j]
      #     b[i] = b[i] - m*b[k];
  
      # for i in range(k+1):
      m = A[k+1:,k]/A[k,k];
      A[k+1:,k:n] = A[k+1:,k:n] - A[k,k:n]*m[:,np.newaxis];
      b[k+1:] = b[k+1:] - m*b[k];

    aux = self.retroSub(A,b,n);

    end= timer()

    return aux, self.erro(A,b,aux),end-start;

  #Decomposição LU
  def decLU(self,A=None,b=None,n=None):
    A,b,n = self.__verifyParameters(A,b,n)
    m = 0
    #Matrizes Extra
    L = np.eye(n,dtype=np.float64());
    U = np.zeros((n,n),dtype=np.float64());

    start = timer()

    # for k in range(0,n):
    #    for i in range(k+1,n):
    #     m = A[i][k]/A[k][k];
    #     L[i][k] = m;
    #     for j in range(k,n):
    #        A[i][j] = A[i][j] - A[k][j]*m;
    # U = A;

    for k in range(0,n):
      m = A[k+1:,k]/A[k,k];
      L[k+1:,k] = m
      A[k+1:] = A[k+1:] - A[k]*m[:,np.newaxis];
      
    U = A;
    subst = self.subst(L,b,n)
    retrosub = self.retroSub(U,subst,n)

    end = timer()

    return retrosub, self.erro(U,b,retrosub), end-start
  
  #Cholesky

  def cholesky(self,A=None,b=None,n=None):
    A,b,n = self.__verifyParameters(A,b,n);

    #Condiçoes: Matriz Positiva
    if(self.__verificaSePositivaDefinida(A) == False):
      return "Matriz Inválida";
   
    G = np.zeros((n,n),dtype = np.float64())

    start =timer()

    # for j in range(0,n):
    #   sum = 0.
    #   for k in range(j):
    #     sum +=G[j][k]**2;
    #   G[j][j] = np.sqrt((A[j][j] - sum));

    #   for i in range(j+1,n):
    #     sum = 0.
    #     for k in range(j):
    #       sum+= G[i][k]*G[j][k]
    #     G[i][j] = (A[i][j]-sum)/G[j][j]

    for j in range(0,n):
      sum = 0.
      
      sum =np.sum(G[j][:j]**2);
      G[j][j] = np.sqrt((A[j][j] - sum));

      sum = 0.
      
      for i in range(j+1,n):
        sum = 0.
        
        sum =np.sum(G[i][:j]*G[j][:j])
        G[i][j] = (A[i][j]-sum)/G[j][j]
      # sum = 0.
      # sum = np.sum(G[:j+1][:j]*G[j][:j])
      # G[:j+1][j] = (A[:j+1][j]-sum)/G[j][j]
        
    Gt = G.transpose()

    subst = self.subst(G,b,n)
    retrosub = self.retroSub(Gt,subst,n)

    end = timer()

    return retrosub, self.erro(Gt,b,retrosub),end-start

  # Faz o produto entre o vetor e a matriz
  def erro(self,A,sol,gerado):
    return np.max(np.abs(np.dot(A,gerado)-sol))
  
  def __verifyParameters(self,A,b,n):
    if((A is None) and (b is None) and (n is None)):
      return self.A.copy(),self.b.copy(),self.n;
    else:
      return A,b,n

  def __verificaSePositivaDefinida(self,A):
      return np.all(np.linalg.eigvals(A) > 0)
  
  # JACOBI

def jacobi(A,b,x,E,n,max):
  #verifica se o chute inicial é nulo
  if x is None:
    x = np.zeros(len(A[0]),dtype=np.float64())
  x_ = np.copy(x)

  start = timer()

  # #vetor dos elementos da diagonal principal de A
  # D = np.diag(A)
  # #matriz com os elementos fora da diagonal principal de A
  # R = A - np.diagflat(D)

  #itera por n vezes
  cont = 0

  while (cont < max):
    
    for i in range (0,n):
      k = np.zeros(max,dtype=np.float64())
      for j in range (0,n):
        if (i != j):
          k[cont] = k[cont] + A[i][j] * x[j]

      x_[i] = (b[i] - k[cont]) / A[i][i]

    erro = np.linalg.norm(x_-x,np.inf)
    if (erro < E):
      end = timer()
      return x_, end-start

    cont += 1
    x = np.copy(x_)
    # x = (b - np.dot(R,x)) / D

  end = timer()  
  return x_, (end-start)

# GAUSS-SEIDEL

def gaussSeid(A,b,x,E,n,max):

  #verifica se o chute inicial é nulo
  if x is None:
    x = np.zeros(len(A[0]),dtype=np.float64())
  x_ = np.copy(x)
  
  start = timer()
  cont=0
  
  while True:
    cont+=1
    for i in range(n):
      
      y = b[i]

      for j in range(0,n):
        if (i!= j):
          y-= A[i][j] * x_[j]
      x_[i] = y/A[i][i]
    erro = np.linalg.norm(x_- x,np.inf)
    if(erro<E):
      end = timer()
      return x_ , (end-start)

    if(cont == max):
      end = timer()
      return x_ , (end-start)
      
    x = np.copy(x_)

  end = timer()
  return x, (end-start)