# Basics of Power Electronics

A basic understanding of power electronics is useful for understanding the operation of the photovoltaic emulator, since the OwnTech TWIST board is a reconfigurable power-conversion platform.

This document presents a concise analysis of the **classical non-synchronous Buck converter**. The photovoltaic emulator itself uses a synchronous Buck-based topology, in which the diode is replaced by an actively controlled semiconductor switch. However, the classical Buck converter provides a simpler introduction to the main operating principles, including switching, duty cycle, inductor behavior, and voltage conversion ratio.

# Buck Converter

A **Buck converter**, also called a step-down converter, is a non-isolated DC–DC converter used to obtain a lower DC output voltage from a higher DC input voltage. In a non-isolated topology, the input and output share a common electrical reference and no transformer is used to provide galvanic isolation.

The classical Buck converter considered here consists of an input DC source $V_{\mathit{in}}$, an electronic switch $S$, a freewheeling diode $D_i$, an inductor $L$, an output capacitor $C$, and a resistive load $R$.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/19b25527-5b48-4687-8a06-848ec98ace8a"
    alt="Classical Buck converter"
    width="550"
  />
</p>

The instantaneous output voltage can be separated into an average DC component and a ripple component:

```math
v_{\mathit{out}}(t)
=
V_{\mathit{out}}
+
\widetilde{v}_{\mathit{out}}(t)
```

where $V_{\mathit{out}}$ is the average output voltage and $\widetilde{v}_{\mathit{out}}$ represents the switching-induced output-voltage ripple.

Power converters operate through the high-frequency switching of semiconductor devices. The switching period $T_s$ and switching frequency $f_s$ are related by

```math
f_s = \frac{1}{T_s}
```

The fraction of each switching period during which the main switch conducts is defined by the duty cycle

```math
D
=
\frac{t_{\mathit{on}}}{T_s}
```

with

```math
0 < D \leq 1
```

# Continuous Conduction Mode

The following analysis assumes operation in **continuous conduction mode (CCM)**, meaning that the inductor current remains greater than zero throughout the entire switching period.

The main electrical quantities considered in the analysis are illustrated below.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/923d6b06-5eed-4fb7-9b37-4b06cfddaf7f"
    alt="Electrical variables of the Buck converter"
    width="550"
  />
</p>

The relevant quantities are:

- $V_{\mathit{in}}$ — input DC voltage
- $v_{DS}$ — MOSFET drain–source voltage
- $i_{DS}$ — MOSFET drain–source current
- $v_{GS}$ — MOSFET gate–source control voltage
- $v_{D_i}$ — diode voltage
- $i_{D_i}$ — diode current
- $v_L$ — inductor voltage
- $i_L$ — inductor current
- $V_{\mathit{out}}$ — average output voltage

For a simplified introductory analysis, the output-voltage ripple is initially neglected. The output voltage is therefore assumed to remain constant over one switching period:

```math
v_{\mathit{out}}(t)
\approx
V_{\mathit{out}}
```

## First Operating Interval: Switch ON

During the first interval,

```math
0 < t \leq D T_s
```

the MOSFET is ON and the diode is reverse-biased.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/429124e3-9b83-45d9-9c0b-32e20f53a4f6"
    alt="Buck converter during the MOSFET ON interval"
    width="550"
  />
</p>

Applying Kirchhoff's voltage law gives the inductor voltage

```math
v_L
=
V_{\mathit{in}}
-
V_{\mathit{out}}
```

Using

```math
v_L
=
L\frac{di_L}{dt}
```

the inductor-current slope is

```math
\frac{di_L}{dt}
=
\frac{
V_{\mathit{in}}
-
V_{\mathit{out}}
}{L}
```

and the current during this interval can be written as

```math
i_L(t)
=
I_m
+
\frac{
V_{\mathit{in}}
-
V_{\mathit{out}}
}{L}t
```

where $I_m$ denotes the minimum inductor current within the switching period.

Since a Buck converter operates with

```math
V_{\mathit{out}} < V_{\mathit{in}}
```

the inductor current increases linearly during this interval. The inductor therefore stores energy in its magnetic field.

For ideal switching devices:

```math
i_{D_i}=0
```

```math
v_{D_i}=-V_{\mathit{in}}
```

```math
i_{DS}=i_L
```

```math
v_{DS}=0
```

## Second Operating Interval: Switch OFF

During the second interval,

```math
D T_s < t \leq T_s
```

the MOSFET is OFF and the diode conducts, providing a path for the inductor current.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/4fee3c02-14c5-4b31-a185-bd774424597f"
    alt="Buck converter during the MOSFET OFF interval"
    width="550"
  />
</p>

The inductor voltage is now

```math
v_L
=
-
V_{\mathit{out}}
```

and therefore

```math
\frac{di_L}{dt}
=
-
\frac{
V_{\mathit{out}}
}{L}
```

The inductor current during this interval is

```math
i_L(t)
=
I_M
-
\frac{
V_{\mathit{out}}
}{L}
\left(
t-DT_s
\right)
```

where $I_M$ denotes the maximum inductor current, reached at the end of the first switching interval.

The negative current slope indicates that the inductor releases part of the energy previously stored in its magnetic field.

For ideal switching devices:

```math
i_{D_i}=i_L
```

```math
v_{D_i}=0
```

```math
i_{DS}=0
```

```math
v_{DS}=V_{\mathit{in}}
```

# Switching Waveforms

The following figure illustrates the main voltage and current waveforms of the classical Buck converter operating in CCM over one switching period.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/b952cd59-a616-414a-89d5-dfc7125ce586"
    alt="Buck converter switching waveforms in continuous conduction mode"
    width="350"
  />
</p>

During the ON interval, the positive inductor voltage causes $i_L$ to increase linearly. During the OFF interval, the negative inductor voltage causes $i_L$ to decrease. Under steady-state CCM operation, these two variations compensate each other over every switching period.

# Static Conversion Ratio

The static conversion ratio of the Buck converter is defined as

```math
G
=
\frac{
V_{\mathit{out}}
}{
V_{\mathit{in}}
}
```

In steady state, the inductor current is periodic:

```math
i_L(t)
=
i_L(t+T_s)
```

Using the inductor relation

```math
v_L
=
L\frac{di_L}{dt}
```

and integrating over one complete switching period gives

```math
\int_{t_0}^{t_0+T_s}
v_L(t)\,dt
=
L
\left[
i_L(t_0+T_s)
-
i_L(t_0)
\right]
=
0
```

Therefore, the average inductor voltage over one switching period is zero. This property is known as **inductor volt-second balance**.

For the Buck converter,

```math
\left(
V_{\mathit{in}}
-
V_{\mathit{out}}
\right)
D T_s
-
V_{\mathit{out}}
\left(
1-D
\right)
T_s
=
0
```

Dividing by $T_s$ and rearranging gives

```math
V_{\mathit{out}}
=
D V_{\mathit{in}}
```

Consequently, the ideal static conversion ratio is

```math
G
=
\frac{
V_{\mathit{out}}
}{
V_{\mathit{in}}
}
=
D
```

with

```math
0 < D \leq 1
```

The ideal Buck converter therefore behaves as a step-down converter: the average output voltage is controlled by the duty cycle and cannot exceed the input voltage.

# Relation to the PV Emulator

The OwnTech TWIST power stage used by the photovoltaic emulator is based on a **synchronous Buck topology** rather than the classical diode-based converter analyzed above.

In a synchronous Buck converter, the freewheeling diode is replaced by an actively controlled semiconductor switch. This reduces conduction losses and enables more flexible power-flow control, while the fundamental relation between duty cycle, inductor voltage, and output-voltage regulation remains similar.

In the photovoltaic emulator, the duty-cycle command $D$ is generated by the embedded voltage controller so that the measured output voltage follows the reference $V_{\mathit{out}}^{\mathit{ref}}$ determined by the emulation algorithm.

For details on how this reference is calculated from the photovoltaic characteristic and the connected load, see **[Strategy.md](Strategy.md)**.
