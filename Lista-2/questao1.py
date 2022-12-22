from sistema_linear import *

#Instancia objeto
n = 10
sist = SistLinear(n)
print("N = " + str(n) + "\n")


#Eliminacação de Gauss(sem pivoteamento)
gauss = sist.eliGauss()
print("Eliminacação de Gauss (sem pivoteamento)")
print("\nVetor solução:")
print(gauss[0])
print("\nErro:")
print(gauss[1])
print("\nTempo de execução:")
print(gauss[2])

#Eliminacação de Gauss(com pivoteamento)
print("\nEliminacação de Gauss (com pivoteamento)")
gauss_piv = sist.pivoteamentoParcial()
print("\nVetor solução:")
print(gauss_piv[0])
print("\nErro:")
print(gauss_piv[1])
print("\nTempo de execução:")
print(gauss_piv[2])

#Decomposição LU
decLU = sist.decLU()
print("\nDecomposição LU")
print("\nVetor solução:")
print(decLU[0])
print("\nErro:")
print(decLU[1])
print("\nTempo de execução:")
print(decLU[2])

#Decomposição de Cholesky
cholesky = sist.cholesky()
print("\nDecomposição de Cholesky")
print("\nVetor solução:")
print(cholesky[0])
print("\nErro:")
print(cholesky[1])
print("\nTempo de execução:")
print(cholesky[2])