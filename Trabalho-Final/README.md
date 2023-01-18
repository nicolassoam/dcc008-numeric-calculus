## ESPECIFICAÇÃO.

---

Implementação do trabalho final proposto para a Disciplina DCC008X - Cálculo Numérico.
---

## OBJETIVOS
Solucionar integrais de um **f(x)** qualquer numericamente através do método de integração de **Gauss-Legendre**.
Para isso, devem ser implementadas o cálculo dos pesos de integração utilizados pelo método para um domínio de integração **[a,b]**.

### ALGORITMOS PROPOSTOS
* **NEWTON RAPHSON**: Método numérico para a solução de sistemas não-lineares. Contém como parâmetros um **chute inicial**, o cálculo de uma matriz **Jacobiana**, formada pelas derivadas parciais de cada função F(wi,ti) e um peso **wi**,**ti**.
* **GAUSS-LEGENDRE**: A partir dos pesos obtidos pelo método de **NEWTON RAPHSON**, realizar o cálculo de integrais de um f(x) qualquer em um domínio de integração arbitário [a,b], dado um número N de pontos para o método.
