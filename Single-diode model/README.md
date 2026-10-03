# Single-Diode Model PV Emulator

This folder contains the implementation of the photovoltaic emulator based on the conventional **single-diode model**.

In this implementation, the photovoltaic $I-V$ characteristic is generated directly from manufacturer-datasheet information. The unknown parameters of the single-diode model are identified numerically, the photovoltaic current is then calculated at selected voltage points, and the resulting characteristic is stored as a piecewise-linear approximation for use by the real-time emulation algorithm.

The general load-identification and operating-point determination procedure used by the emulator is described in **[Strategy.md](../Strategy.md)**.

# Emulator Principle

This implementation enables the OwnTech board to operate as a photovoltaic emulator capable of reproducing the behavior of a PV module under different irradiance and temperature conditions.

To configure the emulator, the user must provide the characteristic quantities normally available in the datasheet of the target module. These quantities are specified at the manufacturer reference conditions, typically **STC** (Standard Test Conditions) or **NOCT** (Nominal Operating Cell Temperature / datasheet reference condition, depending on the module documentation). Once these reference values are defined, the user can specify the irradiance and temperature corresponding to the operating condition to be emulated.

# Single-Diode Photovoltaic Model

The photovoltaic characteristic is generated using the conventional single-diode equivalent model illustrated below.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/f71de800-dbda-4f83-918b-bb4bbc1a45a7"
    alt="Single-diode equivalent model"
    width="500"
  />
</p>

In this model, $I_{\mathit{ph}}$ represents the photogenerated current, $I_s$ the diode saturation current, $R_s$ the series resistance, $R_p$ the parallel resistance, and $A$ the equivalent diode ideality factor of the PV module. The photovoltaic voltage and current are denoted by $V_{\mathit{pv}}$ and $I_{\mathit{pv}}$, respectively.

The photogenerated current is calculated as

```math
I_{\mathit{ph}}
=
I_{\mathit{ph}}^{\mathit{ref}}
\left(
\frac{G}{G^{\mathit{ref}}}
\right)
\left[
1+\alpha
\left(
T-T^{\mathit{ref}}
\right)
\right]
```

where $G$ and $T$ are the irradiance and cell temperature of the operating condition to be emulated, while $G^{\mathit{ref}}$ and $T^{\mathit{ref}}$ correspond to the reference environmental conditions specified in the PV-module datasheet.

The diode saturation current varies with temperature according to

```math
I_s
=
I_s^{\mathit{ref}}
\left(
\frac{T}{T^{\mathit{ref}}}
\right)^3
\exp
\left[
\frac{
n_s q E_G
}{
A k_B
}
\left(
\frac{1}{T^{\mathit{ref}}}
-
\frac{1}{T}
\right)
\right]
```

where $n_s$ is the number of series-connected cells, $q$ is the elementary charge, and $k_B$ is the Boltzmann constant. The semiconductor bandgap energy is calculated using the Varshni equation

```math
E_G
=
E_{G0}
-
\frac{k_1 T^2}{k_2+T}
```

For silicon cells, the values adopted in this implementation are

```math
E_{G0} = 1.166~\mathrm{eV}
```

```math
k_1 = 4.73\times10^{-4}~\mathrm{eV/K}
```

```math
k_2 = 636~\mathrm{K}
```

The resulting photovoltaic-current equation is

```math
I_{\mathit{pv}}
=
I_{\mathit{ph}}
-
I_s
\left[
\exp
\left(
\frac{
q
\left(
R_s I_{\mathit{pv}}+V_{\mathit{pv}}
\right)
}{
A k_B T
}
\right)
-1
\right]
-
\frac{
R_s I_{\mathit{pv}}+V_{\mathit{pv}}
}{
R_p
}
```

Because $I_{\mathit{pv}}$ appears on both sides of the equation, the photovoltaic current cannot be isolated analytically and must therefore be determined numerically.

# Parameter Identification

The single-diode model contains five module-dependent parameters that are not normally provided directly in manufacturer datasheets:

- $I_{\mathit{ph}}^{\mathit{ref}}$ — reference photogenerated current
- $I_s^{\mathit{ref}}$ — reference diode saturation current
- $A$ — equivalent diode ideality factor
- $R_s$ — series resistance
- $R_p$ — parallel resistance

These five parameters are identified from the following nonlinear system:

```math
\begin{cases}

I_{\mathit{ph}}^{\mathit{ref}}
-
I_s^{\mathit{ref}}
\left[
\exp
\left(
\frac{
q R_s I_{\mathit{sc}}^{\mathit{ref}}
}{
A k_B T^{\mathit{ref}}
}
\right)
-1
\right]
-
\frac{
R_s I_{\mathit{sc}}^{\mathit{ref}}
}{
R_p
}
-
I_{\mathit{sc}}^{\mathit{ref}}
=0

\\[1.2em]

I_{\mathit{ph}}^{\mathit{ref}}
-
I_s^{\mathit{ref}}
\left[
\exp
\left(
\frac{
q V_{\mathit{oc}}^{\mathit{ref}}
}{
A k_B T^{\mathit{ref}}
}
\right)
-1
\right]
-
\frac{
V_{\mathit{oc}}^{\mathit{ref}}
}{
R_p
}
=0

\\[1.2em]

I_{\mathit{ph}}^{\mathit{ref}}
-
I_s^{\mathit{ref}}
\left[
\exp
\left(
\frac{
q
\left(
R_s I_{\mathit{mpp}}^{\mathit{ref}}
+
V_{\mathit{mpp}}^{\mathit{ref}}
\right)
}{
A k_B T^{\mathit{ref}}
}
\right)
-1
\right]
-
\frac{
R_s I_{\mathit{mpp}}^{\mathit{ref}}
+
V_{\mathit{mpp}}^{\mathit{ref}}
}{
R_p
}
-
I_{\mathit{mpp}}^{\mathit{ref}}
=0

\\[1.2em]

R_s
+
\frac{
q I_s^{\mathit{ref}} R_p
\left(
R_s-R_p
\right)
}{
A k_B T^{\mathit{ref}}
}
\exp
\left(
\frac{
q R_s I_{\mathit{sc}}^{\mathit{ref}}
}{
A k_B T^{\mathit{ref}}
}
\right)
=0

\\[1.2em]

I_{\mathit{ph}}^{\mathit{ref}}
-
I_s^{\mathit{ref}}
\left\{
\left[
1+
\frac{
q
\left(
V_{\mathit{mpp}}^{\mathit{ref}}
-
R_s I_{\mathit{mpp}}^{\mathit{ref}}
\right)
}{
A k_B T^{\mathit{ref}}
}
\right]
\exp
\left(
\frac{
q
\left(
R_s I_{\mathit{mpp}}^{\mathit{ref}}
+
V_{\mathit{mpp}}^{\mathit{ref}}
\right)
}{
A k_B T^{\mathit{ref}}
}
\right)
-1
\right\}
-
\frac{
2V_{\mathit{mpp}}^{\mathit{ref}}
}{
R_p
}
=0

\end{cases}
```

Because this system is nonlinear, its solution requires an iterative numerical method. In the embedded implementation, the unknown parameters are grouped into the vector

```math
\mathbf{x}
=
\left[
I_{\mathit{ph}}^{\mathit{ref}},
\log_{10}\left(I_s^{\mathit{ref}}\right),
A,
\log_{10}\left(R_s\right),
\log_{10}\left(R_p\right)
\right]^T
```

The logarithmic representation of $I_s^{\mathit{ref}}$, $R_s$, and $R_p$ improves the numerical conditioning of the identification problem while ensuring positive physical values.

The nonlinear system is solved using the **Levenberg–Marquardt algorithm**, with the Jacobian evaluated numerically using forward finite differences.
# $I-V$ Curve Generation

Once the five model parameters have been identified, the photovoltaic characteristic is generated for the irradiance $G$ and temperature $T$ selected by the user.

Since repeatedly solving the full implicit single-diode equation during real-time emulation would increase the computational burden of the controller, the characteristic is computed in advance and stored as a piecewise-linear approximation.

The implementation uses

```math
N_{\mathit{pt}} = 26
```

voltage points with a nonuniform distribution. More points are concentrated around the maximum power point, where the photovoltaic characteristic exhibits stronger nonlinear behavior.

To position this denser region under the selected temperature, the open-circuit voltage is first estimated as

```math
V_{\mathit{oc}}^{\mathit{est}}
=
V_{\mathit{oc}}^{\mathit{ref}}
\left[
1+
\beta
\left(
T-T^{\mathit{ref}}
\right)
\right]
```

The corresponding maximum-power-point voltage is then estimated by

```math
V_{\mathit{mpp}}^{\mathit{est}}
=
V_{\mathit{mpp}}^{\mathit{ref}}
\frac{
V_{\mathit{oc}}^{\mathit{est}}
}{
V_{\mathit{oc}}^{\mathit{ref}}
}
```

For each selected voltage point $V_{\mathit{pv},k}$, the corresponding current $I_{\mathit{pv},k}$ is calculated by solving the implicit single-diode equation using the **Newton–Raphson method**.

Once the voltage–current pairs have been obtained, consecutive points are connected with straight-line segments. For the $k$-th segment,

```math
I_{\mathit{pv},k}^{\mathit{seg}}
\left(
V_{\mathit{pv}}
\right)
=
a_k V_{\mathit{pv}} + b_k
```

with

```math
a_k
=
\frac{
I_{\mathit{pv},k+1}
-
I_{\mathit{pv},k}
}{
V_{\mathit{pv},k+1}
-
V_{\mathit{pv},k}
}
```

and

```math
b_k
=
I_{\mathit{pv},k}
-
a_k V_{\mathit{pv},k}
```

The resulting piecewise-linear characteristic is stored and then used by the real-time emulation algorithm. Consequently, the Levenberg–Marquardt parameter-identification stage and the Newton–Raphson current-calculation stage do not need to be repeated during normal real-time operation.

# Required Inputs

The user does not need to provide the five parameters of the single-diode model directly. They are identified internally from the datasheet quantities and the selected environmental operating conditions.

The following inputs are required:

<div align="center">

<table>
  <thead>
    <tr>
      <th>Parameter</th>
      <th>Description</th>
      <th>Unit</th>
      <th>Source</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>$n_s$</td>
      <td>Number of cells connected in series</td>
      <td>–</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$V_{\mathit{mpp}}^{\mathit{ref}}$</td>
      <td>Maximum-power-point voltage at reference conditions</td>
      <td>V</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$I_{\mathit{mpp}}^{\mathit{ref}}$</td>
      <td>Maximum-power-point current at reference conditions</td>
      <td>A</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$V_{\mathit{oc}}^{\mathit{ref}}$</td>
      <td>Open-circuit voltage at reference conditions</td>
      <td>V</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$I_{\mathit{sc}}^{\mathit{ref}}$</td>
      <td>Short-circuit current at reference conditions</td>
      <td>A</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$G^{\mathit{ref}}$</td>
      <td>Reference irradiance</td>
      <td>W/m²</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$T^{\mathit{ref}}$</td>
      <td>Reference cell temperature</td>
      <td>K</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$\alpha$</td>
      <td>Temperature coefficient of $I_{\mathit{sc}}$</td>
      <td>1/K</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$\beta$</td>
      <td>Temperature coefficient of $V_{\mathit{oc}}$</td>
      <td>1/K</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$G$</td>
      <td>Irradiance to be emulated</td>
      <td>W/m²</td>
      <td>User-defined</td>
    </tr>
    <tr>
      <td>$T$</td>
      <td>Cell temperature to be emulated</td>
      <td>K</td>
      <td>User-defined</td>
    </tr>
  </tbody>
</table>

</div>

The reference quantities correspond to the environmental conditions specified by the manufacturer, typically STC or another datasheet reference condition such as NOCT.

Temperature coefficients must be expressed as fractional values. Therefore, coefficients provided by the manufacturer in `%/K` or `%/°C` must be divided by 100 before being used by the model.

# Test Voltage Used for Load Identification

During the load-identification phase described in **[Strategy.md](../Strategy.md)**, this implementation uses

```math
V_{\mathit{out}}^{\mathit{test}}
=
1.05 \cdot V_{\mathit{oc}}^{\mathit{est}}
```

The factor $1.05$ was selected as an empirical compromise between load-identification accuracy and the voltage applied to the connected device.

> **Safety warning:** This choice must not be interpreted as a general operating requirement. The test voltage and the resulting current must always remain within the admissible limits of the connected device, the OwnTech power-converter platform, and the external DC source. Although the load is initially unknown to the algorithm, the electrical limits of the connected device must be known before the identification procedure is enabled. For sensitive devices, or when these limits cannot be reliably ensured, a lower test voltage or a gradual current-limited identification procedure should be preferred.

# Auxiliary Algorithms

This folder also contains auxiliary MATLAB and Python scripts intended to support the analysis and verification of the single-diode-model implementation.

These scripts are used for tasks such as:

- identifying the five parameters of the single-diode model
- generating and plotting the theoretical $I-V$ characteristic
- building and visualizing the piecewise-linear approximation
- comparing emulation results with theoretical or experimental data

The exact scripts available in this folder can be found in:

- **[Auxiliary Algorithms / MATLAB](./Auxiliary%20Algorithms/MATLAB/)**
- **[Auxiliary Algorithms / Python](./Auxiliary%20Algorithms/Python/)**

# Usage

Before uploading the firmware to the SPIN board, the PV-module datasheet parameters and the desired operating conditions must be configured in `main.cpp`.

The embedded algorithm then:

1. identifies the five single-diode-model parameters using Levenberg–Marquardt
2. generates the photovoltaic characteristic using Newton–Raphson at 26 voltage points
3. constructs and stores the corresponding piecewise-linear approximation
4. executes the real-time emulation strategy described in **[Strategy.md](../Strategy.md)**

For detailed instructions on configuring, compiling, uploading, and operating the emulator, see **[Tutorial.md](../Tutorial.md)**.
