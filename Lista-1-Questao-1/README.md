Seja o problema de valor inicial associado a lei de resfriamento de
Newton: $$\begin{aligned}
\begin{cases}
 \dfrac{d \theta}{d t} = -K(\theta - \theta_m)\\ 
\theta(0) = \theta_0
\end{cases}, 
\quad t\in [0,T]
\end{aligned}$$

A solução exata para o problema [\[N\]](#N){reference-type="eqref"
reference="N"} é dada por: $$
\theta(t) = (\theta_0 - \theta_m) e^{-Kt} + \theta_m$$

-   Apresente as seguintes discretizações por diferenças finitas para
    aproximar o problema [\[N\]](#N){reference-type="eqref"
    reference="N"}:

    -   Euler Explícito:

    -   Euler Implícito;

    -   Crank-Nicolson.

    Pelo método de Euler Explícito, tem-se a seguinte diferenciação:
    $$\begin{aligned}

    \frac{\theta_{n+1}-\theta_{n}}{\Delta t} &=& -K (\theta_n - \theta_m) \Rightarrow\nonumber\\
    \Rightarrow\theta_{n+1}-\theta_n &=& -K\Delta t(\theta_n - \theta_m) \Rightarrow\nonumber\\
    \Rightarrow\theta_{n+1} &=& \theta_n -K\Delta t(\theta_n-\theta_m)
    \end{aligned}$$

    Pelo método de Euler Implícito:
    $$\begin{aligned}
   
    \frac{\theta_{n+1} - \theta_n}{\Delta t} &=& -K(\theta_{n+1}-\theta_m) \Rightarrow \nonumber \\
    \Rightarrow - \theta_n &=& -K\Delta t (\theta_{n+1}-\theta_m) - \theta_{n+1} \Rightarrow \nonumber \\
    \Rightarrow \theta_n &=& K\Delta t (\theta_{n+1}-\theta_m) + \theta_{n+1} \Rightarrow \nonumber \\
    \Rightarrow \theta_n &=& K\Delta t \theta_{n+1} - K\Delta t \theta_m + \theta_{n+1} \Rightarrow \nonumber \\
    \Rightarrow \theta_n &=& \theta_{n+1} (K\Delta t + 1) - K\Delta t \theta_m \Rightarrow \nonumber \\
    \Rightarrow \theta_n + K\Delta t \theta_m &=& \theta_{n+1}(K\Delta t + 1) \Rightarrow \nonumber \\ \nonumber \\
    \Rightarrow \theta_{n+1} &=& \left(\frac{\theta_n + K\Delta t\theta_m}{1+K\Delta t}\right)
    \end{aligned}$$

    Pelo método de Crank-Nicolson:
    $$
    \frac{\theta_{ n+\frac{1}{2} } - \theta_n} {\frac{\Delta t}{2}} = -K(\theta_n - \theta_m)$$

    $$
    \frac{2(\theta_{n+1} - \theta_n)}{\Delta t} = -K (\theta_{n+1} + \theta_n - 2\theta_m)\\$$
