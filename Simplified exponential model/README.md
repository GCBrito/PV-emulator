# Simplified Exponential Model PV Emulator

This folder contains a previous implementation of the photovoltaic emulator based on a **simplified exponential model**.

Unlike the [single-diode model](../Single-diode%20model/), this implementation does not explicitly model the influence of irradiance and temperature on the photovoltaic characteristic. Instead, the $I-V$ curve is generated directly from four characteristic values provided to the algorithm.

The general load-identification and operating-point determination procedure used by the emulator is described in **[Strategy.md](../Strategy.md)**.

# Simplified Exponential Model

The photovoltaic characteristic is represented by the following exponential relation:

```math
I_{\mathit{pv}}
\left(
V_{\mathit{pv}}
\right)
=
I_{\mathit{sc}}^{\mathit{ref}}
\left[
1-
\exp
\left(
\frac{
V_{\mathit{pv}}
-
V_{\mathit{oc}}^{\mathit{ref}}
}{
c
}
\right)
\right]
```

The parameter $c$ is calculated from the maximum-power-point and open-circuit quantities according to

```math
c
=
-
\frac{
V_{\mathit{oc}}^{\mathit{ref}}
-
V_{\mathit{mpp}}^{\mathit{ref}}
}{
\ln
\left(
1-
\frac{
I_{\mathit{mpp}}^{\mathit{ref}}
}{
I_{\mathit{sc}}^{\mathit{ref}}
}
\right)
}
```

where:

- $V_{\mathit{oc}}^{\mathit{ref}}$ — open-circuit voltage
- $I_{\mathit{sc}}^{\mathit{ref}}$ — short-circuit current
- $V_{\mathit{mpp}}^{\mathit{ref}}$ — maximum-power-point voltage
- $I_{\mathit{mpp}}^{\mathit{ref}}$ — maximum-power-point current

These four quantities define the photovoltaic characteristic reproduced by the emulator. They are normally obtained from the manufacturer datasheet at a specified reference condition, such as STC or NOCT.

Because irradiance and temperature are not explicit inputs of the simplified exponential model, the implementation does not internally recalculate the characteristic when these environmental conditions change. To emulate another operating condition, the corresponding values of $V_{\mathit{oc}}$, $I_{\mathit{sc}}$, $V_{\mathit{mpp}}$, and $I_{\mathit{mpp}}$ must therefore be supplied.

# $I-V$ Curve Generation

Once the simplified exponential model has been defined from the four datasheet parameters, the photovoltaic characteristic is calculated in advance and stored as a piecewise-linear approximation for use by the real-time emulation algorithm.

The implemented version uses

```math
N_{\mathit{pt}} = 11
```

predefined voltage points distributed along the $I-V$ characteristic, with a greater concentration of points around the maximum-power and open-circuit regions. For each voltage point $V_{\mathit{pv},k}$, the corresponding current $I_{\mathit{pv},k}$ is calculated directly from the simplified exponential model.

Consecutive voltage–current pairs are then connected by straight-line segments. For the $k$-th segment,

```math
I_{\mathit{pv},k}^{\mathit{seg}}
\left(
V_{\mathit{pv}}
\right)
=
a_k V_{\mathit{pv}} + b_k
```

The resulting piecewise-linear characteristic is stored and subsequently used by the real-time emulation strategy described in **[Strategy.md](../Strategy.md)**.

# Required Inputs

To configure this implementation, the user must provide the following characteristic quantities for the PV module:

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
      <td>$V_{\mathit{mpp}}^{\mathit{ref}}$</td>
      <td>Maximum-power-point voltage at the reference condition</td>
      <td>V</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$I_{\mathit{mpp}}^{\mathit{ref}}$</td>
      <td>Maximum-power-point current at the reference condition</td>
      <td>A</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$V_{\mathit{oc}}^{\mathit{ref}}$</td>
      <td>Open-circuit voltage at the reference condition</td>
      <td>V</td>
      <td>Datasheet</td>
    </tr>
    <tr>
      <td>$I_{\mathit{sc}}^{\mathit{ref}}$</td>
      <td>Short-circuit current at the reference condition</td>
      <td>A</td>
      <td>Datasheet</td>
    </tr>
  </tbody>
</table>

</div>

All four quantities must correspond to the same environmental operating condition.

# Test Voltage Used for Load Identification

During the load-identification phase described in **[Strategy.md](../Strategy.md)**, this implementation uses

```math
V_{\mathit{out}}^{\mathit{test}}
=
1.10 \cdot V_{\mathit{oc}}^{\mathit{ref}}
```

> **Safety warning:** The test voltage and resulting current must always remain within the admissible limits of the connected device, the OwnTech power-converter platform, and the external DC source. The factor $1.10$ should not be interpreted as a general operating requirement for other implementations.

# Auxiliary Algorithm

This folder contains the following MATLAB script:

- `tracer_simplified_exponential_model.m` — generates the simplified exponential $I-V$ characteristic and can be used to visualize the emulated operating points, test points, and load lines.

# Usage

Before uploading the firmware to the SPIN board, the four PV-module parameters must be configured in `main.cpp`.

The embedded algorithm then:

1. calculates the exponential-model parameter $c$
2. generates the photovoltaic current at 11 predefined voltage points
3. constructs the corresponding piecewise-linear approximation
4. executes the real-time emulation strategy described in **[Strategy.md](../Strategy.md)**

For detailed instructions on configuring, compiling, uploading, and operating the emulator, see **[Tutorial.md](../Tutorial.md)**.
