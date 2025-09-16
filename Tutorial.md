# _Test's configuration_

To carry out the laboratory experiments, the setup shown in the following figure was implemented:   

<p align="center">
<img width="500" height="400" alt="Test's configuration" src="https://github.com/user-attachments/assets/f6368711-5379-4dcc-9017-c8cffc4b75b9" />
</p>

It consists of the following elements:

- A DC power supply (_V<sub>in</sub>_).  
- The **OwnTech** board, running the Emulator's algorithm.  
- A variable resistive load.  
- An ammeter to measure the output current (_I<sub>out</sub>_).  
- A voltmeter to measure the output voltage (_V<sub>out</sub>_).

# _First steps_

The first step is to download and install the numerical environment required to use the OwnTech board.  To do so, follow the [**OwnTech Tutorial**](https://docs.owntech.org/latest/core/docs/environment_setup/).

# _Installing_  

Once PlatformIO is installed, download the `main.cpp` and `app.ini` files from either emulator version (the [simplified exponential model](https://github.com/GCBrito/PV-emulator/tree/main/Simplified%20exponential%20model/SPIN%20Firmware) or the [single-diode model](https://github.com/GCBrito/PV-emulator/tree/main/Single-diode%20model/SPIN%20Firmware)) and replace the original `main.cpp` and `app.ini` files in the `core/src` directory.

# _Running_ 

Once `main.cpp` is in place, click **PlatformIO: Build** (bottom status bar) to compile the project. Ensure the build completes with **no errors**.

<p align="center">
<img width="1365" height="99" alt="VS Code (Built)" src="https://github.com/user-attachments/assets/dd96c6d4-b47d-45ed-98d5-03f01d42003a" />
</p>

Then, with the OwnTech board connected to your computer via USB, click **PlatformIO: Upload** to flash the algorithm to the **SPIN** board.

<p align="center">
<img width="1365" height="96" alt="VS Code (Serial Monitor)" src="https://github.com/user-attachments/assets/c19289a4-922b-4e5a-bead-54c964498b89" />
</p>

To operate the OwnTech board, open **PlatformIO: Serial Monitor**.  Press **E** to enter **Emulator Mode**.  Press **I** to stop the emulation (Idle Mode).

<p align="center">
<img width="1365" height="96" alt="VS Code (Serial Monitor)" src="https://github.com/user-attachments/assets/37ac49fb-9c28-4813-939c-0ef9333b5922" />
</p>

# _Safety_

To avoid hazardous currents or voltages, it is important to limit the input voltage and current that the DC source can deliver. When applying these limits, the value of _V<sub>in</sub>_ must remain higher than the _V<sub>OC</sub>_ of the PV module to be emulated.

It is also important to highlight that the **OwnTech** board has a maximum output current limitation of 16 A (considering the two parallel-connected Buck converters). This constraint must be taken into account when selecting the PV module to be emulated: in particular, the short-circuit current (_I<sub>SC</sub>_)must remain below this value.  

Furthermore, the limits of the external sensors (ammeter and voltmeter), as well as those of the connected load, must also be considered to ensure that the peripheral equipment is not subjected to potentially hazardous operating conditions.  
