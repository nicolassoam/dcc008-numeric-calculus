from item_b import *

#QUESTÃO C

print("Erro - Euler Explícito: ", erro_exp(time, thM, K, dt, tt, 0))
print("Erro - Euler Implícito: ", erro_impl(time, thM, K, dt, tt, 0))
print("Erro - Crank-Nicolson:  ", erro_nicol(time, thM, K, dt, tt, 0), '\n')

#Gráficos
def plot_exact(legend):
  plt.plot(tt,exact_sol(th0, thM, K, tt),'--')
  plt.legend([legend,'Exata'])
  title = "N = 25, dt = 2.0833333333333335"
  plt.title(title)
  plt.grid()
  plt.show()
  
t = np.linspace(timeI,time,25)
tt = np.linspace(timeI,time,25)

plt.plot(t,exp_euler(time, thM, K, dt, 25)[0])
plot_exact('Euler Explícito')

plt.plot(t,impl_euler(time, thM, K, dt,25)[0])
plot_exact('Euler Implícito')

plt.plot(t,crank_nicol(time, thM, K, dt,25)[0])
plot_exact('Crank_Nicolson')