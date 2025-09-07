# PV-emulator
This repository provides algorithms and resources to use an OwTech board as a **photovoltaic (PV) emulator**.  
It includes the following folders:

- **Single-diode model** — the most complete version of the PV emulator.  
- **Simplified exponential model** — a simpler version of the PV emulator.  
- **Experimental Data** — all datasets from laboratory tests.  
- **Emulators comparison** — a MATLAB algorithm to evaluate and compare the performance of all emulators.  

# OwnTech

[OwnTech](https://owntech.io) is a company based in Toulouse whose mission is to **democratize power electronics** through open-source technologies.  
The implementation of this PV emulator using the OwnTech solution is motivated by its strong commitment to openness, ensuring compatibility and access to **open-hardware** and **open-software** solutions.  

The CNRS-associated foundation supports the creation of a community where users can share implemented code as well as modifications to the hardware and software designs for specific applications.

## Hardware Overview

The **PCB** (Printed Circuit Board) used for the PV emulator consists of two distinct parts:

- **SPIN (Control Stage)** — handles control and measurement.
- **TWIST (Power Stage)** — provides the power conversion stage.

All electronic schematics related to these boards are available on the foundation’s GitHub repository [OwnTech Foundation GitHub](https://github.com/owntech-foundation)
  
  <img width="2048" height="1365" alt="image" src="https://github.com/user-attachments/assets/7391a637-109c-41bc-a8a8-1d0e5023c9b4" />

## SPIN Board

The **SPIN board** integrates an **STM32 microcontroller** that:  
- Generates the **PWM (Pulse Width Modulation)** signals required to drive the switches.  
- Processes measurements from the sensors embedded in the power stage.  
- Provides real-time access to input and output **current** and **voltage** values.

  <img width="1335" height="717" alt="image" src="https://github.com/user-attachments/assets/e2c93b84-2b84-4ae3-9168-859036231f40" />


## TWIST Board

The **TWIST board** can be configured to operate in three main topologies:  
- Synchronous **Buck converter**  
- Synchronous **Boost converter**  
- **Single-phase inverter**

  <img width="545" height="458" alt="image" src="https://github.com/user-attachments/assets/46f266ec-21d0-4aaf-af63-a85be75b0c3b" />


### Operating Ranges

The TWIST converter supports the following input/output ranges:

- **V<sub>high</sub>**: 12 – 100 V  
- **V<sub>low</sub>**: 12 – 72 V  
- **I<sub>high, low</sub>**: up to 8 A per side  

These ranges define the **power limits of the PV emulator** described in this repository.  



