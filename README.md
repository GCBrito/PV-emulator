# Open-source Photovoltaic Emulator

This repository provides algorithms and resources to use an OwnTech board as a **photovoltaic (PV) emulator**.

It includes the following folders:

- **[Single-diode model](Single-diode%20model/)** — the most recent implementation of the PV emulator, based on the conventional single-diode photovoltaic model.
- **[Simplified exponential model](Simplified%20exponential%20model/)** — a previous implementation of the PV emulator based on a simplified exponential model.
- **[Tests Data](Tests%20Data/)** — datasets collected during laboratory tests.
- **[Emulators comparison](Emulators%20comparison/)** — scripts used to evaluate and compare the performance of the different emulator implementations.
- **[Extras](Extras/)** — articles, lectures, and additional resources related to the PV emulator.

To better understand the algorithms presented in this repository, a set of text files has been organized. The following reading order is recommended for readers who are not yet familiar with this work:

1. **[Main README](README.md)** — introduces the concept of PV emulation and presents the OwnTech platform.
2. **[Strategy](Strategy.md)** — explains the emulation strategy adopted by the proposed PV emulator.
3. **[Single-diode model README](Single-diode%20model/)** — describes the single-diode model for PV modules and its implementation in the emulator.
4. **[Tutorial](Tutorial.md)** — provides instructions on how to configure and operate the emulator.
5. **[Simplified exponential model README](Simplified%20exponential%20model/)** *(optional)* — describes the simplified exponential model and its implementation.
6. **[Basics of Power Electronics](Basics%20of%20Power%20Electronics.md)** *(optional)* — presents a classical Buck-converter analysis for readers who are not yet familiar with power electronics.
7. **[Extras](Extras/)** *(optional)* — contains articles, lectures, and additional material related to the PV emulator.

# PV Emulator

A PV module is a system composed of semiconductor materials capable of converting solar energy into electricity. From an electrical perspective, when environmental conditions are sufficient and a load is connected to a PV module, the voltage and current supplied to the load are determined by the interaction between the photovoltaic characteristic and the connected load. This behavior can be represented in the current–voltage ($I-V$) plane, as illustrated below.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/50484e57-8e17-4bc5-af89-bb46c07745dc"
    alt="I-V characteristic and resistive load line"
    width="650"
  />
</p>

In the $I-V$ plane, the red curve represents the photovoltaic characteristic, described by the photovoltaic voltage $V_{\mathrm{pv}}$ and current $I_{\mathrm{pv}}$, under given irradiance and temperature conditions. The blue line represents the load line associated with a resistive load $R_L$. According to Ohm's law, the load current is given by

```math
I_L(V_{\mathrm{out}})
=
\frac{V_{\mathrm{out}}}{R_L}
```

and the slope of the resistive load line is therefore

```math
\frac{1}{R_L}
```

When a resistive load is connected directly to a PV module, the electrical operating point is determined by the intersection between the photovoltaic $I-V$ characteristic and the load line. In the proposed emulator, the voltage coordinate of this intersection is denoted by

```math
V_{\mathrm{out}}^{\ast}
```

and corresponds to the output-voltage value that the emulator must reproduce for the identified load.

A **photovoltaic emulator (PVE)** is a system designed to reproduce the electrical behavior of a real PV module without requiring a physical photovoltaic panel. For a given photovoltaic characteristic and connected load, the emulator must establish the corresponding operating point. The emulator implemented in this repository is based on the OwnTech platform and can be represented by the following simplified diagram.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/88b33ebe-240c-4d69-b906-b29596ad4287"
    alt="Block diagram of the photovoltaic emulator"
    width="550"
  />
</p>

Unlike a real photovoltaic module, the operating point is not established naturally by the emulator. A dedicated control strategy is therefore required to determine $V_{\mathrm{out}}^{\ast}$ from the connected load and the emulated photovoltaic characteristic and to regulate the converter output accordingly. The general strategy adopted in this project is described in **[Strategy.md](Strategy.md)**. It should also be emphasized that a PV emulator does not convert solar energy into electricity; instead, it relies on an external electrical supply, referred to as the **DC source**, which provides the energy delivered to the connected load.

# OwnTech

[OwnTech](https://owntech.io) develops open-hardware and open-software solutions intended to make power electronics more accessible for education, research, and prototyping. The implementation of the PV emulators available in this repository using the OwnTech platform is motivated by its open architecture, which provides access to both the hardware design and the embedded software and facilitates modification, experimentation, and reproduction of the proposed systems.

The OwnTech platform used for the PV emulator consists of two main boards: the **SPIN board**, responsible for control, computation, and measurement processing, and the **TWIST board**, which provides the power-conversion stage. The corresponding electronic schematics, hardware projects, and supporting resources are publicly available through the **[OwnTech Foundation GitHub](https://github.com/owntech-foundation)**.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/7391a637-109c-41bc-a8a8-1d0e5023c9b4"
    alt="OwnTech SPIN and TWIST boards"
    width="750"
  />
</p>

## SPIN Board

The **SPIN board** integrates an **STM32G474RE microcontroller** and provides the embedded control and measurement-processing functions required by the photovoltaic emulator. It acquires the voltage and current measurements provided by the power stage, executes the photovoltaic-emulation and control algorithms, and generates the PWM signals used to drive the switches of the TWIST board. USB-C connectivity also enables communication with a host computer for configuration and monitoring, while the real-time emulation and control calculations are executed directly on the embedded microcontroller.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/090d9a0e-13fb-443f-be09-a5b10a50856c"
    alt="OwnTech SPIN control board"
    width="450"
  />
</p>

## TWIST Board

The **TWIST board** provides the power-conversion stage of the system. It is a reconfigurable bidirectional power converter rated at **300 W**, with two low-side channels operating from **12 V to 72 V** and rated at **8 A each**, and one high-side channel operating from **10 V to 110 V** and rated at **16 A**. Integrated voltage and current sensors provide the measurements required by the embedded control system.

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/46f266ec-21d0-4aaf-af63-a85be75b0c3b"
    alt="OwnTech TWIST power board"
    width="450"
  />
</p>

The TWIST board can be configured for different power-converter topologies depending on the application. In the photovoltaic-emulator implementations provided in this repository, it is used as a Buck-based power-conversion stage. The exact topology, configuration, and control strategy associated with each emulator are described in the corresponding implementation documentation.

# Reproducibility

This repository is intended to make the photovoltaic-emulator implementations, experimental data, numerical tools, and associated documentation publicly available. Implementation-specific source code and explanations are provided in the **[Single-diode model](Single-diode%20model/)** and **[Simplified exponential model](Simplified%20exponential%20model/)** folders, while experimental datasets and comparison tools are provided in **[Tests Data](Tests%20Data/)** and **[Emulators comparison](Emulators%20comparison/)**.

The electronic schematics and open-hardware design files of the OwnTech platform are maintained separately by OwnTech and are publicly available through the **[OwnTech Foundation GitHub](https://github.com/owntech-foundation)**. A versioned release of this repository will be used to identify the exact version associated with the corresponding scientific publication.
