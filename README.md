
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
  
  <img width="1500" height="800" alt="image" src="https://github.com/user-attachments/assets/7391a637-109c-41bc-a8a8-1d0e5023c9b4" />

### SPIN Board

The **SPIN board** integrates an **STM32 microcontroller** that:  
- Generates the **PWM (Pulse Width Modulation)** signals required to drive the switches.  
- Processes measurements from the sensors embedded in the power stage.  
- Provides real-time access to input and output **current** and **voltage** values.

  <img width="1335" height="717" alt="image" src="https://github.com/user-attachments/assets/e2c93b84-2b84-4ae3-9168-859036231f40" />


### TWIST Board

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

# PV Emulator

A Photovoltaic (PV) module is a system composed of semiconductor materials capable of converting solar energy into electricity. From an electrical perspective, when environmental conditions are sufficient, connecting a load to a PV panel automatically subjects it to a DC (Direct Current) voltage and current. This principle is best illustrated by examining the current–voltage (I–V) plane, as shown in the figure below:

<img width="3248" height="1948" alt="I-V plan" src="https://github.com/user-attachments/assets/93737631-7bcd-4a1c-89a6-baca75cddae9" />

In this plane, the red curve represents the **I–V characteristic**, which describes the electrical behavior of a PV module (each module has its own curve) under given temperature and irradiance conditions. The blue line, known as the **load line**, corresponds to the resistive load R connected to the module. Its slope, defined by Ohm’s Law, is equal to the inverse of the resistance value. When a resistive load is connected to a PV panel, the operating point of the system is determined by the intersection between the module’s I–V curve and the load line. This point, denoted as (_V*_,_I*_), specifies the voltage _V*_ across the load and the current _I*_ flowing through it. In a real PV panel, this operating point is reached naturally, without external intervention. 

It is also possible, however, to design a system capable of reproducing this same electrical behavior: a **_PV emulator_**. The emulator proposed in this repository is based on OwnTech technology and can be represented by the following diagram:

<img width="1374" height="483" alt="Emulator" src="https://github.com/user-attachments/assets/88b33ebe-240c-4d69-b906-b29596ad4287" />

A PV emulator is a system that replicates the electrical characteristics of a real solar panel. In other words, for a given resistive load R, it delivers the same voltage and current (_V*_,_I*_) that the load would receive if it were directly connected to a PV module. Unlike a real PV panel, though, the operating point (_V*_,_I*_) cannot be achieved naturally and requires a dedicated control strategy, described in the section *Emulation Strategy*. It should also be emphasized that, unlike an actual PV module, a PV emulator does not convert solar energy into electricity. Instead, it relies on an external electrical supply (referred to as the DC source in the figure) as its power input.






