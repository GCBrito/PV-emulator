This emulator is based on the single-diode mathematical model. To configure the emulation, the user must specify the photovoltaic (PV) module to replicate and define the desired operating point by providing the following input parameters:

- ns: number of cells in series
- np: number of parallel strings
- VMPmod: voltage at the maximum power point
- IMPmod: current at the maximum power point
- VOCmod: open-circuit voltage
- ISCmod: short-circuit current
- Tref: reference temperature
- Sref: reference irradiance
- muICC: temperature coefficient of the short-circuit current
- T: actual temperature
- S: actual irradiance

This folder contains two MATLAB scripts:

- "tracer_single_diode_emulator.m": Determines the five parameters required for the single-diode model by solving a nonlinear system of equations. Additionally, it can plot both the emulated operating points and the PV emulator’s testing points, and compare them to the theoretical I–V curve.
- "comparison_single_diode_emulator_real_pv.m": Compares the emulated I–V curve to experimental data from a real PV module.

*Usage instructions:*

To use the emulator, begin by running "tracer_single_diode_emulator.m" to compute the five model parameters. Once obtained, these parameters must be manually inserted into main.cpp before uploading the firmware to the SPIN board
