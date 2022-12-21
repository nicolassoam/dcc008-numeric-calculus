from item_c import *


#QUESTÃO D

print("delta t variando (letra b)")

print("Erro - Euler Explícito: ", erro_exp(time, thM, K, dt, tt, 50))
print("Erro - Euler Implícito: ", erro_impl(time, thM, K, dt, tt, 50))
print("Erro - Crank-Nicolson:  ", erro_nicol(time, thM, K, dt, tt, 50), '\n')

print("delta t único (letra c) => N = 25")

print("Erro - Euler Explícito: ", erro_exp(time, thM, K, dt, tt, 25))
print("Erro - Euler Implícito: ", erro_impl(time, thM, K, dt, tt, 25))
print("Erro - Crank-Nicolson:  ", erro_nicol(time, thM, K, dt, tt, 25), '\n')
