from funcoes import *


#QUESTÃO B

# N para Euler Implícito
N_euler = [3, 5, 15, 25,50]
# N para Crank-Nicolson
N_crank =  [2, 3, 6, 15,50]

#Gráficos
for n in range(0,5):

  t = np.linspace(timeI,time,N_euler[n])

  #time
  tt = np.linspace(timeI,time,N_euler[n])

  exp_eu = exp_euler(time, thM, K, dt,N_euler[n])
  impl_eu = impl_euler(time, thM, K, dt,N_euler[n])
  plt.plot(t,exp_eu[0],'*')
  plt.plot(t,impl_eu[0],'*')
  plt.plot(tt,exact_sol(th0, thM, K, tt))
  plt.legend(['Euler Explícito','Euler Implícito','Exata'])
  title = "N = " + str(N_euler[n]) + " , dt = " + str(impl_eu[1])
  plt.title(title)
  plt.grid()
  plt.show()

for n in range(0,5):

  t = np.linspace(timeI,time,N_crank[n])

  #time
  tt = np.linspace(timeI,time,N_crank[n])

  #vetor de solucao (32 bits)
  sol = np.zeros(N_crank[n],dtype=np.float32())
  sol[0] = th0
  crank = crank_nicol(time, thM, K, dt,N_crank[n])
  plt.plot(t,crank[0],'*')
  plt.plot(tt,exact_sol(th0, thM, K, tt))
  plt.legend(['Crank-Nicolson','Exata'])
  title = "N = " + str(N_crank[n]) + " , dt = " + str(crank[1])
  plt.title(title)
  plt.grid()
  plt.show()
