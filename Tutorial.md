<img width="1365" height="629" alt="image" src="https://github.com/user-attachments/assets/6bcd7cc4-54b7-435f-9ba6-5dec2697cbe6" /># _Test's configuration_

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

After completing the previous step, copy the `main.cpp` file from your chosen emulator version — either the [**Simplified Exponential Model**](https://github.com/GCBrito/PV-emulator/tree/main/Simplified%20exponential%20model/SPIN%20Firmware) or the [**Single-Diode Model**](https://github.com/GCBrito/PV-emulator/tree/main/Single-diode%20model/SPIN%20Firmware). You can do this by clicking on the **"Copy raw file"** option in GitHub.

<img width="1365" height="628" alt="image" src="https://github.com/user-attachments/assets/27199407-4c56-4c73-9f1e-044866413654" />

Once this is done, paste the copied algorithm into the `main.cpp` file that was generated in the PlatformIO environment after completing OwnTech’s **“Environment Setup”** tutorial.

# _Operating_ 

Once `main.cpp` is in place, click **PlatformIO: Build** (bottom status bar) to compile the project. Ensure the build completes with **no errors**.

<p align="center">
<img width="1365" height="99" alt="VS Code (Built)" src="https://github.com/user-attachments/assets/dd96c6d4-b47d-45ed-98d5-03f01d42003a" />
</p>

Then, with the OwnTech board connected to your computer via USB, click **PlatformIO: Upload** to upload the algorithm to the **SPIN** board.

<p align="center">
<img width="1365" height="96" alt="VS Code (Serial Monitor)" src="https://github.com/user-attachments/assets/c19289a4-922b-4e5a-bead-54c964498b89" />
</p>

To operate the OwnTech board, open **PlatformIO: Serial Monitor**. Press **P** to enter **Power Mode**. Press **E** to enter **Emulator Mode**. Press **I** to stop the emulation (**Idle Mode**).

<p align="center">
<img width="1365" height="96" alt="VS Code (Serial Monitor)" src="https://github.com/user-attachments/assets/37ac49fb-9c28-4813-939c-0ef9333b5922" />
</p>

# _Safety_

To avoid hazardous currents or voltages, it is essential to limit the input voltage and current that the DC source can deliver by properly configuring the source itself.  
It is also important to note that when applying these limits, the value of _V<sub>in</sub>_ must remain higher than the _V<sub>OC</sub>_ of the PV module being emulated, in order to ensure proper operation of the **Emulator Mode**.

In addition, the **OwnTech** board has a maximum output current limitation of **16 A**, considering the two parallel-connected Buck converters.  
This constraint must be taken into account when selecting the PV module to be emulated — in particular, the short-circuit current (_I<sub>SC</sub>_) must remain below this value.

Finally, the limits of the external sensors (ammeter and voltmeter), as well as those of the connected load, must also be respected to ensure that the peripheral equipment is not exposed to potentially hazardous operating conditions.
