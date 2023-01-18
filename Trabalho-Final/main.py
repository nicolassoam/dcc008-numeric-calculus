from methods import *

# DEFINIÇÃO DE PARAMETROS

# Número de pontos
N = 4

# Limites
a,b = -1,1


#Número MÁXIMO de iterações
max = 500
# Pertubação Epsilon
eps = 10**(-8)

# Alpha para regularização da matriz
alpha = 10**(-4)


# CHAMADA DE FUNÇÕES

# Chute Inicial
# x = chuteInicial(a,b,N)

# Integral
# integral = calculaIntegrais(a,b,N)
# print(integral)


# print(x[0])

# Matriz Jacobiana
# jacob = jacobiana(x[0],eps,N,a,b,integral,alpha)
# print(jacob)


# solucao = newtonRaphson(a,b,x[0],-eps,N,max,alpha)
# print(solucao)

integracao = GaussLegendre(a,b,N,eps,alpha,max)
print("Solução da integral de f(x) = sqrt(a**2 + b) aproximada pela quadratura de Gauss Legendre:")
print(integracao)
print("Solução exata:")
exata = integralExata2(a,b)
print(exata)
print("")
print("Erro: ")
print(erro(exata,integracao))
print("")

print("Número de pontos = " + str(N))
print("Limites: ")
print("a = " + str(a) + " b = "+ str(b))