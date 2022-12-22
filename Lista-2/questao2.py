from sistema_linear import *

# EQ DIFERENCIAL RESOLVIDA PELO MÉTODO DAS DIFERENÇAS FINITAS

def equation(n,metodo,is_iterative,E):

  print("n = " + str(n) + " | E = " + str(E) + "\n");

  k = int(np.sqrt(n))
  A = np.zeros((n,n))
  b = 10*np.ones(n)
  u = np.zeros(n)

  for i in range(0,n):
    A[i,i] = -4
  for i in range(0,n-1):
    A[i+1,i] = 1
    A[i,i+1] = 1
    if((i+1)%k ==0):
      A[i+1,i] = 0
      A[i,i+1] = 0
  for i in range(k,n):
    A[i-k,i] = 1
    A[i,i-k] = 1

  A = - A*(k-1)**2

  print("\nMatriz Gerada por método de Diferenças Finitas:\n")
  print(A)

  u = np.linalg.solve(A,b)

  # print("Solução dada por np.solve()")
  # print(u)

  if(is_iterative):
    aux = metodo(A,b,None,E,n,300)
    print("\nSolução dada por métodos iterativos:")
    print(aux[0])
    print("\nTempo de Execução:")
    print(aux[1])
  else:
    aux = metodo(A,b,n)
    print("\nSolução dada por métodos diretos:")
    print(aux[0])
    print("\nTempo de Execução:")
    print(aux[2])
   
  # x = np.linspace(0,1,k)
  # y = np.linspace(0,1,k)

  # X,Y = np.meshgrid(x,y)

  # U = np.zeros((k,k))
  # L = np.zeros((k,k))

  # for i in range(0,k):
  #   for j in range(0,k):
  #     U[i,j] = u[i+j*k]
  #     L[i,j] = aux[0][i+j*k]

  # fig = plt.figure()
  # fig2 = plt.figure()

  # # ax = fig.add_subplot(1,1,1, projection='3d')
  # bx = fig2.add_subplot(1,1,1, projection='3d')

  # # ax.plot_surface(X, Y, U,label='NP.SOLVE()')
  # bx.plot_surface(X, Y, L,label='MÉTODOS IMPLEMENTADOS')

  # plt.show()
 
  return aux

sist = SistLinear()



# n = 81, 289, 1089, 4225, 16641
n = 1089

# tolerância
E = 10**(-5)

print("Jacobi")
M1 = equation(n,jacobi,True,E)
print("\nGauss-Seidel")
M2 = equation(n,gaussSeid,True,E)
print("\nCholesky")
M3 = equation(n,sist.cholesky,False,E)
print("\nDecomposição LU")
M4 = equation(n,sist.decLU,False,E)
print("\nGauss (sem pivoteamento)")
M5 = equation(n,sist.eliGauss,False,E)
print("\nGauss (com pivoteamento)")
M6 = equation(n,sist.pivoteamentoParcial,False,E)