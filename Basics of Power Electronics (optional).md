To better understand how the PV emulator operates, a brief overview on **power electronics** is helpful — especially since the OwnTech board is a power-electronics platform. This section presents a concise analysis of the **classical (non-synchronous) Buck converter** — rather than the synchronous Buck used on the OwnTech board — to cover the core concepts; the operating principles are essentially the same in both cases.

# _Buck Converter_ 

A **Buck (step-down) converter** is a **non-isolated DC–DC** circuit that reduces a higher DC input voltage to a lower DC output voltage. “Non-isolated” means **input and output share the same ground**—there is no transformer providing electrical isolation.
The following figure shows a Buck converter composed of a DC input voltage _V<sub>in</sub>_, an inductor _L_, an electronic switch _S_ (in this overview, a MOSFET was considered), a diode _D<sub>i</sub>_, a capacitor _C_, and a load _R_. The voltage _v<sub>o</sub>_ is the output voltage associated with the load. Note that _v<sub>o</sub>_ consists of two components: a DC (average) component _V<sub>o</sub>_ and an AC component $$\tilde{v}_o
$$, commonly referred to as the **output voltage ripple**.

<p align="center">
<img width="500" height="400" alt="Figura 6" src="https://github.com/user-attachments/assets/19b25527-5b48-4687-8a06-848ec98ace8a" />
</p>

In this types of systems, operation is cyclical and relies on high-frequency **switching** of semiconductor devices (transistors) that act as electronic switches. The **switching period** _T<sub>s</sub>_ and **switching frequency** _f<sub>s</sub>_ describe the repetitive ON/OFF cycle and are related by _f<sub>s</sub> = 1/T<sub>s</sub>_.

## Analysis in CCM 

To understand the general characteristics and behavior of the topology in **continuous conduction mode (CCM)** — meaning that, in steady state, the inductor current never falls to zero — a deeper analysis is required. This analysis describes the voltage and current waveforms in the main components and derives the **static conversion ratio** \(G\), defined as the ratio of the average output to the average input voltage:

$$
G = \frac{\{V}_o}{{V}_{in}}
$$

To conduct this analysis, the variables of interest associated with each element of the examined configuration are first defined : 

<p align="center">
<img width="500" height="400" alt="Figura 7" src="https://github.com/user-attachments/assets/923d6b06-5eed-4fb7-9b37-4b06cfddaf7f" />
</p>

In Figure 10, the variables highlighted in **red** represent:

- _V<sub>in</sub>_ - average input DC voltage
- _v<sub>DS</sub>_ - MOSFET drain–source voltage  
- _i<sub>DS</sub>_ - MOSFET drain–source current  
- _v<sub>GS</sub>_ - MOSFET gate–source (control) voltage
- _v<sub>Di</sub>_ - diode voltage  
- _i<sub>Di</sub>_ - diode current  
- _v<sub>L</sub>_ - inductor voltage  
- _i<sub>L</sub>_ - inductor current  
- _V<sub>o</sub>_ - average output voltage across the load _R_ (output voltage)

> **Modeling note.** For the sake of a simpler analysis, in a first moment, it's assumed that the output voltage is ripple-free, i.e., it remains constant and equal to its average value _V<sub>o</sub>_. This simplification is standard practice in introductory power-electronics analyses.

### First operating interval (MOSFET ON)

The first operating stage of the Boost converter occurs over $0 < t \le &alpha; \ T_s$, where &alpha; is called the duty cycle and is defined as

$$
&alpha; = \frac{t_{\mathrm{on}}}{T_s},
$$

i.e., the fraction of the total switching period $T_s$ during which the MOSFET $S$ conducts. In other words, $t_{\mathrm{on}}$ is the time the switch is ON within one period of the gate–source control signal $v_{GS}$.

During this stage, the **MOSFET conducts** while the **diode is reverse-biased (OFF)** :

<p align="center">
<img width="500" height="400" alt="Figura 8" src="https://github.com/user-attachments/assets/429124e3-9b83-45d9-9c0b-32e20f53a4f6" />
</p>

During this half-cycle, since the diode behaves as an open circuit and the capacitor voltage equals the output voltage, the inductor voltage is:

$$
v_L = V_{in} - V_o.
$$

From $v_L = L \dfrac{di_L}{dt}$, the inductor current evolution in this interval follows:

$$
i_L(t) = \frac{V_{in} - V_o}{L}  t + I_m
$$

From this equation, it follows that the inductor current increases linearly during the first interval because, for a Buck converter, $V_o < V_{in}$ and thus $\tfrac{di_L}{dt} = \tfrac{V_{in}-V_o}{L} > 0$. Physically, the inductor is storing energy in its magnetic field during this stage. In this expression, $I_{m}$ denotes the minimum inductor current within one switching period (i.e., typically $i_L(0)=I_{m}$).

Regarding the commutation elements in this interval (MOSFET ON, diode OFF), it follows that:
- $i_D = 0$ (the diode is reverse-biased),
- $v_D = -V_{in}$,
- $i_{DS} = i_L$ (the MOSFET drain–source current equals the inductor current),
- $v_{DS} = 0$ (an ideal switch in conduction has no voltage drop).

### Second operating interval (MOSFET OFF)

During the time interval $\alpha T_s < t \le T_s$, the converter enters the second operating phase, in which the **diode is ON** and the **MOSFET is OFF**:

<p align="center">
<img width="500" height="400" alt="Figura 9" src="https://github.com/user-attachments/assets/4fee3c02-14c5-4b31-a185-bd774424597f" />
</p>

From the analysis of this circuit, it is possible to determine that the inductor voltage is:

$$
v_L = -V_o.
$$

This leads to the inductor current:

$$
i_L(t) = -\frac{V_o}{L}\big(t - \alpha T_s\big) + I_M
$$

From this equation, it is possible to notice that the inductor current $i_L$ decreases linearly with slope $-V_o/L$. Physically, this means that the inductor releases the energy stored in its magnetic field during the first stage of operation. In this expression, $I_{M}$ denotes the maximum inductor current within one switching period (i.e., typically $i_L(\alpha T_s)=I_{M}$).

For the remaining circuit elements in this interval (diode ON, MOSFET OFF), one obtains:

Regarding the commutation elements in this interval (MOSFET ON, diode OFF), it follows that:
- $i_D = i_L$,
- $v_D = 0$,
- $i_{DS} = 0$,
- $v_{DS} = V_{o}$.

### Waveforms 

The following figure depicts the time evolution of voltages and currents in the Buck topology in CCM over one switching period _T<sub>s</sub>_:

<p align="center">
<img width="250" height="700" alt="waveforms" src="https://github.com/user-attachments/assets/b952cd59-a616-414a-89d5-dfc7125ce586" />
</p>

### Static Conversion Ratio (Buck)

After reviewing the time-domain behavior, an important quantity — the **static conversion ratio** (voltage gain) — can be computed. To proceed, note that in steady state (CCM) the inductor current is **periodic**:

$$
i_L(t) = i_L(t+T_s),\quad \forall t.
$$

Using $v_L = L \dfrac{di_L}{dt}$ and integrating over one period:

$$
\int_{t_0}^{t_0+T_s} v_L(t)dt
= L\int_{t_0}^{t_0+T_s} \frac{di_L}{dt}dt
= L[i_L(t_0+T_s) - i_L(t_0)] = 0
$$

which implies that the **average inductor voltage is zero** over a switching period. Thus, the algebraic sum of the areas A<sup>+</sup> and A<sup>-</sup> illustrated on the previous image over one switching period is zero:

$$
(V_{in}-V_o) \alpha T_s + (-V_o)(T_s - \alpha T_s) = 0.
$$

Dividing by T<sub>s</sub> and rearranging:

$$
(V_{in}-V_o) \alpha - V_o (1-\alpha) = 0
\Rightarrow
V_o = \alpha V_{in}.
$$

Therefore, the static gain is

$$
G = \frac{V_o}{V_{in}} = \alpha , \qquad 0 < \alpha \le 1.
$$

This expresses the step-down nature of the Buck converter: the output voltage is equal to or lower than the input.


