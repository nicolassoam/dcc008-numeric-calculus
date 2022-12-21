# Lista 1 de Cálculo Numérico - Parte 1
# Maria Luísa Riolino Guimarães (202165563C)
# Nicolas Soares Martins (202165559C)

import numpy as np
import matplotlib.pyplot as plt
import math

#solução exata.
def exact_sol(th0, thM, K, tt):
  return (th0 - thM)*np.exp((-K)*tt) + thM

#calcula delta t
def calc_dt(time, timeI):
   return (time - timeI)/(time-1)

#Euler Explícito
def exp_euler(time, thM, K, dt, N):

  #se houver refinamento
  if(N > 0):
    dt = (time-timeI)/(N-1)
    #vetor de solucao (64 bits)
    sol = np.zeros(N,dtype=np.float32())
    sol[0] = th0
    for n in range(0,(N-1)):
      sol[n+1] = dt*(-K*(sol[n] - thM)) + sol[n]
  else:
    #vetor de solucao (64 bits)
    sol = np.zeros((time),dtype=np.float32())
    sol[0] = th0
    for n in range(0,(time-1)):
        sol[n+1] = dt*(-K*(sol[n] - thM)) + sol[n]
  return sol,dt

#Euler Implícito.
def impl_euler(time, thM, K, dt, N):
  if(N > 0):
    dt = (time-timeI)/(N-1)
    sol = np.zeros(N,dtype=np.float32())
    sol[0] = th0
    for n in range(0,(N-1)):
      sol[n+1] = (sol[n] + K*dt*thM)/(1+K*dt)
  else:
    sol = np.zeros((time),dtype=np.float32())
    sol[0] = th0
    for n in range(0,(time-1)):
      sol[n+1] = (sol[n] + K*dt*thM)/(1+K*dt)
  return sol,dt 

#Crank-Nicolson
def crank_nicol(time, thM, K, dt,N):
  if(N > 0):
    dt = (time-timeI)/(N-1)
    sol = np.zeros(N,dtype=np.float32())
    sol[0] = th0 
    for n in range(0,(N-1)):
      sol[n+1] = ((-K*dt/2)*(sol[n]-2*thM)+sol[n])/((K*dt)/2 + 1)
  else:
    sol = np.zeros((time),dtype=np.float32())
    sol[0] = th0
    for n in range(0,(time-1)):
        sol[n+1] = ((-K*dt/2)*(sol[n]-2*thM)+sol[n])/((K*dt)/2 + 1)
  return sol,dt  

#Cálculo do Erro (Euler Explícito)
def erro_exp(time, thM, K, dt, tt, N):
  if(N>0):
    tt = np.linspace(timeI,time,N)
  return np.max(np.abs(exact_sol(th0, thM, K, tt) - exp_euler(time, thM, K, dt, N)[0]))

#Cálculo do Erro (Euler Implícito)
def erro_impl(time, thM, K, dt, tt, N):
  if(N>0):
    tt = np.linspace(timeI,time,N)
  return np.max(np.abs(exact_sol(th0, thM, K, tt) - impl_euler(time, thM, K, dt, N)[0]))

#Cálculo do Erro (Crank-Nicolson)
def erro_nicol(time, thM, K, dt, tt, N):
   if(N>0):
    tt = np.linspace(timeI,time,N)
   return np.max(np.abs(exact_sol(th0, thM, K, tt) - crank_nicol(time, thM, K, dt, N)[0]))



#constante K
K = 0.035871952

#condicao inicial
th0 = 99.

#theta m
thM = 27.

# Intervalo
time = 50
timeI = 0

dt = calc_dt(time, timeI)