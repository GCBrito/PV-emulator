This emulator is based on a simplified exponential mathematical model. To configure the emulation, the user must specify the photovoltaic (PV) module to be replicated and define the desired operating point by providing the following input parameters:

- Vmp: voltage at the maximum power point
- Imp: current at the maximum power point
- Voc: open-circuit voltage
- Isc: short-circuit current

This folder contains one MATLAB script:

- "tracer_simplified_exponential_model.m": Plots the emulated operating points and the PV emulator’s testing points, and compares them to the theoretical I–V curve.

*Usage instructions:*

To use the emulator, the user must manually insert the parameters listed above into main.cpp before uploading the firmware to the SPIN board.