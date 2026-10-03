# Emulation Strategy

The physical configuration of the photovoltaic emulator is illustrated in the following schematic. The DC source provides the electrical energy required by the emulator, while the OwnTech power converter regulates the electrical quantities applied to the connected load. The source-side voltage and current are denoted by $V_{\mathit{in}}$ and $I_{\mathit{in}}$, respectively, while the output voltage and current are denoted by $V_{\mathit{out}}$ and $I_{\mathit{out}}$.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/f4c03241-04e5-4d4b-b164-15f3085446b1" 
    alt="Physical configuration of the photovoltaic emulator"
    width="550"
  />
</p>

The emulation strategy relies on the photovoltaic current–voltage characteristic and on the electrical behavior imposed by the connected load. Its objective is to determine the operating point that would naturally result if the same load were connected to the photovoltaic module being emulated and then regulate the output of the power converter so that this operating condition is reproduced.

# Photovoltaic $I-V$ Characteristic

The electrical behavior of a photovoltaic module is commonly represented by its current–voltage characteristic, in which the photovoltaic voltage and current are denoted by $V_{\mathit{pv}}$ and $I_{\mathit{pv}}$, respectively. A typical $I-V$ characteristic is illustrated below.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/325602d3-1189-4b32-a344-fd1664ef48dd"
    alt="Photovoltaic I-V characteristic"
    width="550"
  />
</p>

Four characteristic quantities are commonly provided in photovoltaic-module datasheets: the open-circuit voltage $V_{\mathit{oc}}$, the short-circuit current $I_{\mathit{sc}}$, the maximum-power-point voltage $V_{\mathit{mpp}}$, and the maximum-power-point current $I_{\mathit{mpp}}$. The maximum power point (MPP) corresponds to the operating condition at which the product of photovoltaic voltage and current is maximized. The $I-V$ characteristic depends on the photovoltaic module and on its operating conditions, particularly irradiance and cell temperature.

This repository contains two different implementations for generating the photovoltaic characteristic:

- **[Simplified exponential model](Simplified%20exponential%20model/)**
- **[Single-diode model](Single-diode%20model/)**

The mathematical formulation and numerical procedures used to generate the photovoltaic characteristic are specific to each model and are therefore described in their corresponding documentation. Nevertheless, both implementations ultimately provide a numerical representation of the photovoltaic $I-V$ characteristic composed of linear segments. The emulation strategy described below operates on this piecewise-linear representation and is therefore common to both implementations.

# Power-Converter Configuration

The photovoltaic emulator is implemented using the OwnTech SPIN control board and TWIST power stage. In the configuration adopted for the emulator, the two low-side channels of the TWIST board are connected in parallel and operated as a two-phase interleaved synchronous Buck converter.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/375fe16a-b458-409b-a89d-ef0d6607a71c"
    alt="Interleaved synchronous Buck converter topology implemented on the TWIST board"
    width="650"
  />
</p>

The SPIN board executes the embedded control algorithm and generates the duty-cycle command $D$ applied to the power converter. An inner voltage-control loop regulates the output voltage, while the outer emulation algorithm identifies the connected load and determines the corresponding operating point on the photovoltaic characteristic.

The two low-side channels are measured independently by the OwnTech sensors. The instantaneous measured output voltage is calculated as the arithmetic mean of the two channel-voltage measurements

```math
V_{\mathit{out}}^{\mathit{OT}}
=
\frac{
V_{\mathit{out},1}^{\mathit{OT}}
+
V_{\mathit{out},2}^{\mathit{OT}}
}{2}
```

while the measured total output current is obtained by summing the two channel-current measurements

```math
I_{\mathit{out}}^{\mathit{OT}}
=
I_{\mathit{out},1}^{\mathit{OT}}
+
I_{\mathit{out},2}^{\mathit{OT}}
```

The superscript $OT$ denotes quantities obtained from the sensors integrated into the OwnTech platform.

# Control Architecture

The general control architecture of the photovoltaic emulator is illustrated below.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/59c4b581-5b8e-4870-b7a3-4f4934727ce5"
    alt="Block diagram of the closed-loop control system"
    width="750"
  />
</p>

The DC source supplies the interleaved Buck converter, whose output is connected to the load. The OwnTech sensors provide the measured output voltage $V_{\mathit{out}}^{\mathit{OT}}$ and current $I_{\mathit{out}}^{\mathit{OT}}$. These quantities are processed by the outer emulation algorithm to identify the connected load and determine the corresponding photovoltaic operating point.

The resulting voltage reference $V_{\mathit{out}}^{\mathit{ref}}$ is compared with the measured output voltage. The corresponding error is processed by the inner voltage controller, which generates the duty-cycle command $D$ applied to the PWM modulator.

# PV-Emulation Sequence

When the PV-emulation mode is activated, the electrical characteristics of the connected load are initially unknown to the emulation algorithm. The system therefore begins with a **load-identification phase**, during which the output-voltage reference is initialized to a predefined test voltage denoted by $V_{\mathit{out}}^{\mathit{test}}$. The exact definition of this voltage depends on the selected emulator implementation and is provided in the corresponding model documentation.

Applying $V_{\mathit{out}}^{\mathit{test}}$ to the connected load produces a voltage–current pair that can be used to identify its electrical behavior. Since the instantaneous sensor measurements may contain switching ripple and noise, the emulation algorithm does not directly use individual samples. Instead, the measured output voltage and current are accumulated over consecutive averaging windows before being used by the outer emulation algorithm.

The real-time control task acquires measurements with a sampling period

```math
T_c = 100~\mu\mathit{s}
```

while the measurements are accumulated over averaging windows of approximately

```math
T_w = 500~\mathit{ms}
```

corresponding to approximately

```math
N_w = 5000
```

samples per averaging window.

For the $w$-th averaging window, the averaged output voltage is defined as

```math
\overline{V_{\mathit{out}}^{\mathit{OT}}}[w]
=
\frac{1}{N_w}
\sum_{r=wN_w}^{(w+1)N_w-1}
V_{\mathit{out}}^{\mathit{OT}}[r]
```

and the averaged output current is

```math
\overline{I_{\mathit{out}}^{\mathit{OT}}}[w]
=
\frac{1}{N_w}
\sum_{r=wN_w}^{(w+1)N_w-1}
I_{\mathit{out}}^{\mathit{OT}}[r]
```

where $r$ denotes the sampling instant of the real-time control loop and $w$ denotes the averaging-window index.

# Quasi-Steady-State Validation

After the test voltage has been applied, the output may require some time to reach a sufficiently stable condition. Therefore, the load-identification procedure is completed only after the averaged output current satisfies a quasi-steady-state criterion.

The relative variation between two consecutive averaged current values is calculated as

```math
er_I
=
\frac{
\left|
\overline{I_{\mathit{out}}^{\mathit{OT}}}[w]
-
\overline{I_{\mathit{out}}^{\mathit{OT}}}[w-1]
\right|
}{
\left|
\overline{I_{\mathit{out}}^{\mathit{OT}}}[w-1]
\right|
}
```

and the measurements are considered sufficiently stable when

```math
er_I < \varepsilon_I
```

with

```math
\varepsilon_I = 0.10
```

This criterion prevents the load-identification procedure from being completed while the output current is still undergoing a significant transient.

# Load Identification

Once the quasi-steady-state condition is satisfied, the most recent averaged output voltage and current are used to estimate the equivalent resistive load.

The load resistance is given by

```math
R_L
=
\frac{
\overline{V_{\mathit{out}}^{\mathit{OT}}}
}{
\overline{I_{\mathit{out}}^{\mathit{OT}}}
}
```

The corresponding load line in the $I-V$ plane is therefore expressed as

```math
I_L(V_{\mathit{out}})
=
\frac{V_{\mathit{out}}}{R_L}
```

This line passes through the origin and represents the electrical behavior of the connected resistive load. The load-identification procedure allows the emulator to determine $R_L$ directly from measured electrical quantities, so the resistance of the connected load does not need to be provided to the emulator in advance.

# Operating-Point Determination

Once $R_L$ has been identified, the corresponding load line is intersected with the piecewise-linear photovoltaic characteristic generated by the selected model.

For the $k$-th linear segment, the photovoltaic current can be represented by

```math
I_{\mathit{pv},k}^{\mathit{seg}}(V)
=
a_k V + b_k
```

where $a_k$ and $b_k$ denote the slope and intercept of the segment.

The desired operating point must simultaneously belong to the photovoltaic characteristic and to the load line. Therefore, its voltage coordinate $V_{\mathit{out}}^{\ast}$ must satisfy

```math
I_{\mathit{pv},k}^{\mathit{seg}}
\left(
V_{\mathit{out}}^{\ast}
\right)
=
I_L
\left(
V_{\mathit{out}}^{\ast}
\right)
```

Substituting the two line equations gives the candidate intersection voltage

```math
V_{\mathit{out}}^{\ast}
=
-
\frac{
b_k
}{
a_k-\frac{1}{R_L}
}
```

The intersection is considered valid only when $V_{\mathit{out}}^{\ast}$ lies within the voltage interval associated with the corresponding segment. Once a valid intersection is found, $V_{\mathit{out}}^{\ast}$ becomes the desired output voltage associated with the operating point of the emulated photovoltaic module.

The voltage reference generated by the outer emulation algorithm can therefore be summarized as

```math
V_{\mathit{out}}^{\mathit{ref}}
=
\begin{cases}
V_{\mathit{out}}^{\mathit{test}}
&
\text{during load identification}
\\
V_{\mathit{out}}^{\ast}
&
\text{during PV emulation}
\end{cases}
```

# Voltage Control

The voltage reference $V_{\mathit{out}}^{\mathit{ref}}$ generated by the outer emulation algorithm is tracked by an inner voltage-control loop. The measured output voltage $V_{\mathit{out}}^{\mathit{OT}}$ is compared with the reference, and the resulting error is processed by the controller to generate the duty-cycle command $D$.

The voltage-control error is expressed as

```math
e_V[r]
=
V_{\mathit{out}}^{\mathit{ref}}[r]
-
V_{\mathit{out}}^{\mathit{OT}}[r]
```

During the load-identification phase, the controller regulates the output voltage toward $V_{\mathit{out}}^{\mathit{test}}$. Once the photovoltaic operating point has been calculated, the reference is updated to $V_{\mathit{out}}^{\ast}$ and the converter is regulated toward this new value.

During normal PV emulation, $V_{\mathit{out}}^{\mathit{ref}}$ remains equal to $V_{\mathit{out}}^{\ast}$ until a change in the connected load is detected.

# Load-Change Detection

After the initial operating point has been established, the emulator continues monitoring the connected load using the averaged output-voltage and output-current measurements.

A candidate load resistance is calculated according to

```math
R_L^{\mathit{cand}}
=
\frac{
\overline{V_{\mathit{out}}^{\mathit{OT}}}
}{
\overline{I_{\mathit{out}}^{\mathit{OT}}}
}
```

and compared with the resistance associated with the previously accepted operating point, denoted by $R_L^{\mathit{prev}}$.

The corresponding relative variation is

```math
er_{R_L}
=
\frac{
\left|
R_L^{\mathit{cand}}
-
R_L^{\mathit{prev}}
\right|
}{
R_L^{\mathit{prev}}
}
```

A possible load change is detected when

```math
er_{R_L} > \varepsilon_R
```

with

```math
\varepsilon_R = 0.02
```

To prevent isolated disturbances or measurement noise from being interpreted as an actual load change, the condition must remain satisfied for

```math
N_{\mathit{det}} = 2
```

consecutive evaluations.

If the condition is not satisfied persistently, the detection counter is reset and normal PV emulation continues. If the condition is confirmed, the algorithm returns to the load-identification phase by restoring the voltage reference to $V_{\mathit{out}}^{\mathit{test}}$. The load is then identified again, and a new operating point $V_{\mathit{out}}^{\ast}$ is calculated.

# Emulation Sequence in the $I-V$ Plane

The complete load-identification and operating-point update sequence is illustrated below.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/405709a5-dd1e-472e-b49a-4b787c7181d2"
    alt="Load-line-based operating-point determination and update sequence"
    width="700"
  />
</p>

The sequence can be described by considering an initially connected resistive load $R_{L1}$ followed by a change to a second load $R_{L2}$:

1. The voltage reference is initialized to $V_{\mathit{out}}^{\mathit{test}}$. The resulting averaged output voltage and current define the test point from which $R_{L1}$ is identified.

2. The load line corresponding to $R_{L1}$ is intersected with the photovoltaic $I-V$ characteristic, and the voltage reference is updated to

```math
V_{\mathit{out}}^{\ast}(R_{L1})
```

3. If the connected load changes from $R_{L1}$ to $R_{L2}$, the voltage reference initially remains equal to $V_{\mathit{out}}^{\ast}(R_{L1})$, while the output current changes according to the new electrical condition.

4. Once the load-change criterion is satisfied, the emulator returns to the load-identification phase by applying $V_{\mathit{out}}^{\mathit{test}}$.

5. The new load $R_{L2}$ is identified, its load line is intersected with the photovoltaic characteristic, and the new operating voltage is calculated

```math
V_{\mathit{out}}^{\ast}(R_{L2})
```

The emulator then resumes normal operation at the updated photovoltaic operating point.

# Implementation-Specific Parameters

The strategy described above is common to both photovoltaic-emulator implementations available in this repository. However, some numerical procedures and parameters depend on the selected photovoltaic model.

These implementation-specific details include:

- the mathematical model used to generate the photovoltaic characteristic
- the number and distribution of points used to discretize the characteristic
- the procedure used to calculate the photovoltaic current
- the definition of $V_{\mathit{out}}^{\mathit{test}}$
- model-specific initialization procedures and numerical parameters

Further details are provided in the corresponding documentation:

- **[Simplified exponential model](Simplified%20exponential%20model/)**
- **[Single-diode model](Single-diode%20model/)**

For practical instructions on configuring and operating the emulator, see **[Tutorial.md](Tutorial.md)**.
