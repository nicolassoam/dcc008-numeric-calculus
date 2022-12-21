from funcoes import *

eps = 0.1;

result = TDMA(a,b,c,d,49,eps,h);
tt = np.linspace(0,1,50);

#print(result);
plt.plot(tt,result,'*-');
plt.plot(t,exact_solution(t,eps,h));
plt.legend(['Aproximada','Exata']);
plt.title("Solução exata x aproximada (ε = 0.1)")
plt.grid();
plt.show();

eps = 0.01;

result = TDMA(a,b,c,d,49,eps,h);
tt = np.linspace(0,1,50);

#print(result);
plt.plot(tt,result,'*-');
plt.plot(t,exact_solution(t,eps,h));
plt.legend(['Aproximada','Exata']);
plt.title("Solução exata x aproximada (ε = 0.01)")
plt.grid();
plt.show();

eps = 0.001;

result = TDMA(a,b,c,d,49,eps,h);
tt = np.linspace(0,1,50);

#print(result);
plt.plot(tt,result,'*-');
plt.plot(t,exact_solution(t,eps,h));
plt.legend(['Aproximada','Exata']);
plt.title("Solução exata x aproximada (ε = 0.001)")
plt.grid();
plt.show();

eps = 0.0001;

result = TDMA(a,b,c,d,49,eps,h);
tt = np.linspace(0,1,50);

#print(result);
plt.plot(tt,result,'*-');
plt.plot(t,exact_solution(t,eps,h));
plt.legend(['Aproximada','Exata']);
plt.title("Solução exata x aproximada (ε = 0.0001)")
plt.grid();
plt.show();