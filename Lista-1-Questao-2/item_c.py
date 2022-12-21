from funcoes import *

base = 0.1;

resultados = np.empty((4,5))

erros = np.zeros(5);
ref = np.zeros(5);

# Cálculo do erro para cada ɛ
for j in range(0,4):
  eps = base * pow(10, -j);
  
  # Cálculo em cada malha
  for i in range(1,6):
    n_el = pow(4,i);
    h = 1/n_el;
    print(n_el)
    #diagonal inferior
    a = np.zeros(n_el-1, dtype=np.float64());

    #diagonal central
    b = np.zeros(n_el, dtype=np.float64());

    #diagonal superior
    c = np.zeros(n_el-1, dtype=np.float64());

    #vetor solucao
    d = np.zeros(n_el, dtype=np.float64());

    t = np.linspace(0,1,n_el);
    tt = np.linspace(0,1,n_el);

    sol = exact_solution(t,eps,h);
    # erros[i-1] = np.log10(np.max(np.abs(sol - TDMA(a,b,c,d,n_el-1,eps,h))));
    # erros[i-1] = np.log10(np.linalg.norm(sol - TDMA(a,b,c,d,n_el-1,eps,h)));
    # ref[i-1] = -np.log10(h);
    erros[i-1] = np.max(np.abs(sol - TDMA(a,b,c,d,n_el-1,eps,h)));
    ref[i-1] = -np.log(h)
    print(erros[i-1])

  resultados[j] = erros;
  nn = np.logspace(0,1,5);
  plt.title("ɛ = " + str(eps))
  # plt.xscale("log",nonposx="clip")
  # plt.yscale("log",nonposy="clip")
  plt.plot(-np.log(ref),np.log(erros));
  plt.grid();
  plt.show();

nn = np.logspace(0,1,5,base=10);
plt.title("Comparativo do erro entre solução exata e aproximada")
# plt.xscale("log",nonposx="clip")
# plt.yscale("log",nonposy="clip")
plt.plot(ref,np.log(resultados[0]),'*-');
plt.plot(ref,np.log(resultados[1]),'^-');
plt.plot(ref,np.log(resultados[2]));
plt.plot(ref,np.log(resultados[3]));

plt.legend(['ɛ = 0.1', 'ɛ = 0.01', 'ɛ = 0.001', 'ɛ = 0.0001']);
plt.grid();
plt.show();