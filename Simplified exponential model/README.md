# _Emulator's Principle_

This repository presents a C++ algorithm that allows the OwnTech board to perfom as PV emulator that can reproduce the behavior of any PV module, but only at the operating points **STC** (Standard Test Conditions — 1000 W/m² solar irradiance, AM1.5 spectrum and cell temperature of 25 °C) or **NOTC** (Normal Operating Test Conditions — 800 W/m² solar irradiance, AM1.5 spectrum, cell temperature of 20 °C and wind speed of 1 m/s). This limitation arises because the emulator is based on the **simplified exponential model**, a mathematical model that reconstructs the I–V characteristic curve of a PV module using only four parameters provided in the manufacturer’s datasheet: _V<sub>OC</sub>_, _I<sub>SC</sub>_, _V<sub>MPP</sub>_, and _I<sub>MPP</sub>_. However, since datasheets typically list these values only at STC or NOTC, the emulator cannot accurately represent the electrical behavior of a PV module under different irradiance or temperature conditions.

# _Simplified exponencial model_

The previously mentioned mathematical model can be described by the following two equations:

$$
I(V) = I_{SC}\left(1 - e^{\tfrac{V - V_{OC}}{c}}\right)
$$

$$
c = -\frac{V_{OC} - V_{MPP}}{\ln\left(1 - \tfrac{I_{MPP}}{I_{SC}}\right)}
$$

Where:  

- _V<sub>OC</sub>_ — open-circuit voltage [V]  
- _I<sub>SC</sub>_ — short-circuit current a [A]  
- _V<sub>MPP</sub>_ — voltage at the maximum power point [V]  
- _I<sub>MPP</sub>_ — current at the maximum power point [A]

Through these equations, the I–V characteristic curve of any PV module can be determined. The simplified exponential model differs from the [single-diode model](https://github.com/GCBrito/PV-emulator/tree/main/Single-diode%20model) because it does not account for the influence of temperature and irradiance on the PV panel. As a result, the emulator can only operate at STC or NOTC conditions, which represents a limited range of operation.

# _Algorithms_

The algorithms in this repository enable the implementation of the simplified exponencial model on the PV emulator.  To use this emulator, the user must specify the PV module to be replicated and define the desired operating point by providing the following input parameters, which are always available in manufacturers’ datasheets:

- Vmp: voltage at the maximum power point [V] 
- Imp: current at the maximum power point [A] 
- Voc: open-circuit voltage [V] 
- Isc: short-circuit current [A] 

This folder contains one **MATLAB** script:

- "tracer_simplified_exponential_model.m": Plots the emulated operating points as well as the PV emulator’s testing points and the load lines, and compares them to the theoretical I–V curve.

# _Usage instructions_

**To use the emulator**, the user must manually insert the parameters listed above into main.cpp before uploading the firmware to the SPIN board. 
