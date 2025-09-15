
This repository provides algorithms and resources to use an OwnTech board as a **photovoltaic (PV) emulator**.  It includes the following folders:

- **Simplified exponential model** — a simplified version of the PV emulator.
- **Single-diode model** — the most complete version of the PV emulator.   
- **Experimental Data** — datasets collected from laboratory tests.  
- **Emulators comparison** — a MATLAB algorithm to evaluate and compare the performance of the different emulators.  

To better understand the algorithms presented in this repository, a set of text files has been organized.  
The following reading order is recommended for readers who are not yet familiar with this work:

1) [**Main README**](https://github.com/GCBrito/PV-emulator/tree/main) — introduces the concept of a PV emulator and presents the OwnTech technology.  
2) [**Strategy**](https://github.com/GCBrito/PV-emulator/blob/main/Strategy.md) — describes the proposed emulation strategy.  
3) [**Simplified exponential model README**](https://github.com/GCBrito/PV-emulator/tree/main/Simplified%20exponential%20model) — explains the simplified exponential model and how the corresponding emulator works.  
4) [**Single-diode model README**](https://github.com/GCBrito/PV-emulator/tree/main/Single-diode%20model) — explains the single-diode model and how the corresponding emulator works.  
5) [**Tutorial**](https://github.com/GCBrito/PV-emulator/blob/main/Tutorial.md) — provides instructions on how to use the emulator in a laboratory environment.
 
# _PV Emulator_

A Photovoltaic (PV) module is a system composed of semiconductor materials capable of converting solar energy into electricity. From an electrical perspective, when environmental conditions are sufficient, connecting a load to a PV panel automatically subjects it to a DC (Direct Current) voltage and current. This principle is best illustrated by examining the current–voltage (I–V) plane, as shown in the figure below:

<p align="center">
<img width="500" height="1000" alt="I-V plan" src="https://github.com/user-attachments/assets/93737631-7bcd-4a1c-89a6-baca75cddae9" />
</p>
  
In this plane, the red curve represents the **I–V characteristic**, which describes the electrical behavior of a PV module (each module has its own curve) under given temperature and irradiance conditions. The blue line, known as the **load line**, corresponds to the resistive load R connected to the module. Its slope, defined by Ohm’s Law, is equal to the inverse of the resistance value. When a resistive load is connected to a PV panel, the operating point of the system is determined by the intersection between the module’s I–V curve and the load line. This point, denoted as (_V*_,_I*_), specifies the voltage _V*_ across the load and the current _I*_ flowing through it. In a real PV panel, this operating point is reached naturally, without external intervention. 

It is also possible, however, to design a system capable of reproducing this same electrical behavior: a **PV emulator**. The emulator proposed in this repository is based on OwnTech technology and can be represented by the following diagram:

<p align="center">
<img width="500" height="400" alt="Emulator" src="https://github.com/user-attachments/assets/88b33ebe-240c-4d69-b906-b29596ad4287" />
</p>

A PV emulator is a system that replicates the electrical characteristics of a real solar panel. In other words, for a given resistive load R, it delivers the same voltage and current (_V*_,_I*_) that the load would receive if it were directly connected to a PV module. Unlike a real PV panel, though, the operating point (_V*_,_I*_) cannot be achieved naturally and requires a dedicated control strategy, described in [Strategy](https://github.com/GCBrito/PV-emulator/blob/main/Strategy.md). It should also be emphasized that, unlike an actual PV module, a PV emulator does not convert solar energy into electricity. Instead, it relies on an external electrical supply (referred to as the DC source in the figure) as its power input.

# _OwnTech_

[OwnTech](https://owntech.io) is a company based in Toulouse whose mission is to **democratize power electronics** through open-source technologies.  
The implementation of this PV emulator using the OwnTech solution is motivated by its strong commitment to openness, ensuring compatibility and access to **open-hardware** and **open-software** solutions.  

The CNRS-associated foundation supports the creation of a community where users can share implemented code as well as modifications to the hardware and software designs for specific applications.

## Hardware Overview

The **PCB** (Printed Circuit Board) used for the PV emulator consists of two distinct parts:

- **SPIN (Control Stage)** — handles control and measurement.
- **TWIST (Power Stage)** — provides the power conversion stage.

All electronic schematics related to these boards are available on the foundation’s GitHub repository [OwnTech Foundation GitHub](https://github.com/owntech-foundation)

<p align="center">
<img width="1000" height="500" alt="image" src="https://github.com/user-attachments/assets/7391a637-109c-41bc-a8a8-1d0e5023c9b4" />
</p>

### SPIN Board

The **SPIN board** integrates an STM32 microcontroller that generates the PWM (Pulse Width Modulation) signals required to drive the switches present in the TWIST board, processes measurements from the sensors embedded in the Power Stage, and provides real-time access to input and output current and voltage values.

<p align="center">
<img width="400" height="300" alt="Adobe Express - file" src="https://github.com/user-attachments/assets/090d9a0e-13fb-443f-be09-a5b10a50856c" />
</p>

### TWIST Board

The **TWIST board** can be configured to operate in three main topologies:  
- Synchronous **Buck converter**  
- Synchronous **Boost converter**  
- **Single-phase inverter**

<p align="center">
<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/46f266ec-21d0-4aaf-af63-a85be75b0c3b" />
</p>

The TWIST converter supports the following input/output ranges:

- **V<sub>high</sub>**: 12 – 100 V  
- **V<sub>low</sub>**: 12 – 72 V  
- **I<sub>high, low</sub>**: up to 8 A per side  

These ranges define the **power limits of the PV emulator** described in this repository.  






