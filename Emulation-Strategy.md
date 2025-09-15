The physical configuration of the PV emulator described here is illustrated in the following schematic, where **_V<sub>in</sub>_** denotes the input voltage imposed by the DC source, **_I<sub>in</sub>_** the input current drawn from the DC source, **_V<sub>s</sub>_** the output voltage imposed by the emulator on the load, and **_I<sub>s</sub>_** the output current supplied by the emulator to this same load:

<p align="center">
<img width="500" height="400" alt="Emulator_parameters" src="https://github.com/user-attachments/assets/601a5e40-206e-4b3b-bcdc-aef52cd53dc9" />
</p>

The PV emulator proposed in this repository relies heavily on the understanding of the **current–voltage (I–V) characteristic curve** of photovoltaic (PV) modules. For this reason, it is important to examine this curve in more detail.

## _I–V Characteristic Curve_  

As previously explained, the electrical behavior of a photovoltaic (PV) panel is generally represented by its I-V curve, which can be illustrated in the next figure:

<p align="center">
<img width="800" height="500" alt="Courbe caractéristique" src="https://github.com/user-attachments/assets/67c87c0e-fbee-4427-a6e2-9b1c867d8de6" />
</p>

In this figure, four key parameters can be identified:

- **_V<sub>oc</sub>_** - Open-circuit voltage  
- **_I<sub>sc</sub>_** - Short-circuit current  
- **_V<sub>MPP</sub>_** - Voltage at the maximum power point  
- **_I<sub>MPP</sub>_** - Current at the maximum power point  

These parameters mainly depend on the intrinsic properties of the PV cells. However, external factors such as **shading, temperature, and environmental conditions** can also affect the overall performance of the module. The **maximum power point (MPP)**, represented by the coordinates (_V<sub>MPP</sub>_, _I<sub>MPP</sub>_) in the I–V plane, corresponds to the operating point at which the **product of current and voltage is maximized**. 

This curve can be modeled using mathematical functions, and several models have been proposed in the literature. In this context, this repository considers two main mathematical models for I–V curves:

- [**Simplified exponential model**](https://github.com/GCBrito/PV-emulator/tree/main/Simplified%20exponential%20model)  
- [**Single-diode model**](https://github.com/GCBrito/PV-emulator/tree/main/Single-diode%20model)

Regardless of which of these two models is considered, both are capable of reproducing the complete I–V characteristic curve of a specific PV module using only a few parameters that are readily available in the datasheet.  

In this context, and with the goal of using the I–V curve in a numerical environment (for the emulation process), an algorithm was developed to reproduce the complete I–V characteristic curve by combining one of the previously mentioned mathematical models with linear interpolation.

## _Emulation Strategy_  

As presented before, the PV emulator proposed in this repository is based on **OwnTech technology** and is specifically configured to operate as two parallel-connected synchronous Buck converters (more information about this configuration can be found in [OwnTech's GitHub](https://github.com/owntech-foundation/examples/blob/main/TWIST/DC_DC/buck_voltage_mode/README.md)):

<p align="center">
<img width="600" height="500" alt="Board" src="https://github.com/user-attachments/assets/03e188e7-65ce-43b2-bf82-c58b0233bb4c" />
</p>

The operation of the emulator relies on controlling the duty cycle **α**, defined as the ratio between the conduction time _t<sub>on</sub>_ of the electronic switches and the total switching period _T<sub>s</sub>_. The duty cycle is managed by a C++ algorithm pre-programmed in the **SPIN** control board, enabling the TWIST to operate as a PV emulator.

### Algorithm  

The algorithm was developed in a **PlatformIO** environment within **Visual Studio Code** and then uploaded to the SPIN board via USB connection. It runs in loop and provides three distinct operating modes:

- **Power Mode (key 'P')** - The board operates as two conventional parallel synchronous Buck converters regulated by a PID controller.  
  - Open-loop mode: the PID tracks a predefined duty cycle reference.  
  - Closed-loop mode: the PID regulates the output to follow a voltage reference.

- **Idle Mode (key 'I')** - Conversion is disabled, and no power is delivered to the load.  

- **Emulator Mode (key 'E')** - The board reproduces the electrical behavior of a PV module supplying a resistive load R. In this mode, the user must first provide the characteristic parameters available in the datasheet of the target PV module. From these values, the algorithm locally approximates the I–V curve through linear interpolation.  

The detailed operation of the **Emulator Mode** is illustrated by the following flowchart, which highlights the main tasks performed by the algorithm. It is important to note that the execution period of these tasks is fixed at **500 µs**.

<p align="center">
<img width="300" height="800" alt="Algorigrame" src="https://github.com/user-attachments/assets/71ec4446-cb0c-4c0c-ab82-129461845094" />
</p>

The nomenclature used in the flowchart is presented below:

- _**V<sub>s0</sub> , I<sub>s0</sub>**_ - Initial voltage and current test-point, used for calculating the load line.
- _**V<sub>s</sub> , I<sub>s</sub>**_ - Output voltage and current, measured across the load by the TWIST sensors.    
- _**R**_ - Load resistance, calculated from V<sub>s</sub> and I<sub>s</sub>.  
- _**V<sub>s*</sub>**_ , _**I<sub>s*</sub>**_ - Voltage and current corresponding to the operating point on the I–V characteristic of the emulated PV module.
- _**α**_ - Duty cycle, control variable used for tuning the output voltage.  

As illustrated by the flowchart, once the user activates the **Emulator Mode** — and assuming a load _R<sub>1</sub>_ connected to the emulator terminals — the duty cycle _α_ is adjusted so that the output voltage _V<sub>s</sub>_ reaches the reference value _V<sub>s0</sub>_, establishing an operating point referred to as the **test point**.  The goal is to used the sensors of the TWIST board are used to measure the voltage _V<sub>s0</sub>_ and the current _I<sub>s0</sub>_ imposed on the load at the test point.

To reduce oscillations caused by noise and ripple, the algorithm records the sensor data from the board for a duration of 500 µs, and then computes their time-averaged values. This process provides more reliable estimates of _V<sub>s</sub>_ and _I<sub>s</sub>_.  

After two successive averaged values of the current _I<sub>s</sub>_ are obtained (i.e., after two full iterations of the previous step), the program enters a **waiting loop**. This loop continues until the relative error between the two most recent _I<sub>s</sub>_ measurements is less than **2%**. This condition ensures that steady-state operation has been reached.  

Once this condition is satisfied, the most recent measured values of _V<sub>s</sub>_ and _I<sub>s</sub>_ (ideally equal to _V<sub>s0</sub>_ and _I<sub>s0</sub>_) are used to calculate the equivalent load resistance _R<sub>1</sub>_, which is unknown a priori. Graphically, this resistance corresponds to the **slope of the load** line passing through the origin of the I–V plane and the test point measured when _V<sub>s</sub> = V<sub>s0</sub>_.  

The mathematical expression of this load line is:

$$
I_{s,1} = \frac{1}{R_{1}} \cdot V_{s,1}
$$

When the emulator mode is activated with a load **R<sub>1</sub>** connected, the duty cycle is initialized to a predefined value **α<sub>0</sub>**.  
To reduce oscillations caused by noise and ripple, the algorithm collects sensor data from the TWIST board for **500 µs** and computes time-averaged values of **V<sub>s</sub>** and **I<sub>s</sub>**.  

Once two successive averaged values of **I<sub>s</sub>** are obtained (after two full iterations), the program enters a waiting loop until the absolute difference between the two values is less than **0.45 A**, ensuring steady-state operation.  

At this point, the last measured values of **V<sub>s</sub>** and **I<sub>s</sub>** are used to compute the equivalent load resistance **R<sub>1</sub>**:  

\[
I_{s,1} = \frac{1}{R_1} \cdot V_{s,1}
\]  

From the load line equation, the algorithm identifies the intersection point (_V<sub>s,1</sub><sup>*</sup>_, _I<sub>s,1</sub><sup>*</sup>_) between this line and the linearly-approximated I–V characteristic of the studied PV module. This point corresponds to the voltage and current that would be imposed across the load _R<sub>1</sub>_ if it were connected to a real PV module.  

Once this point is determined, the voltage reference is updated, and the PID controller adjusts the duty cycle until _V<sub>s</sub>_ reaches the new target value _V<sub>s,1</sub><sup>*</sup>_ (closed loop).

If, at a later stage, the load _R<sub>1</sub>_ connected to the system is replaced by a new load _R<sub>2</sub>_, the algorithm automatically adapts to this change.  This adaptability is ensured by a periodic task executed every 500 µs, during which a candidate resistance, _R<sub>candidate</sub>_, is estimated from the most recent measurements of _V<sub>s</sub>_ and _I<sub>s</sub>_. The relative error between _R<sub>candidate</sub>_ and the reference resistance _R<sub>1</sub><sup>*</sup>_ is then computed according to the following formula :

$$
Error = \frac{|R_{2}-R_{1}^{\ast}|}{R_{1}^{\ast}}
$$

If the error calculated according to this equation exceeds **2%** for two consecutive cycles of the periodic task (i.e., a total of 1000 µs), the algorithm detects a load change and reinitializes the voltage reference to duty cycle to _V<sub>s0</sub>_, thereby **restarting** the adaptation process. It is important to note that this persistent error condition over two cycles prevents transient variations caused by noise or ripple from being mistaken for an actual load change.

To facilitate the understanding of the **Emulator Mode** logic, the following figure illustrates the operation of the emulator in the I–V plane:

<p align="center">
<img width="500" height="1000" alt="Emulator_Working" src="https://github.com/user-attachments/assets/a981126c-f72d-4fa1-8829-520670ffce35" />
</p>

In this image, the numbered points represent the key steps of the **Emulator Mode** operation:

1. The initial voltage reference _V<sub>s0</sub>_ is applied on the load, generating a load line passing through the origin.  
2. The PID then adjusts the duty cycle (_α<sub>1</sub><sup>*</sup>_) to reach the operating point (_V<sub>s,1</sub><sup>*</sup>_, _I<sub>s,1</sub><sup>*</sup>_) on the I–V curve.  
3. When the load changes, the operating point shifts to point (3). Since _α<sub>1</sub><sup>*</sup>_ remains constant, the reference voltage _V<sub>s,1</sub><sup>*</sup>_ does not change; only the current varies in response to the new load.  
4. If the relative error between the previous load and the candidate load exceeds 2% for two consecutive cycles, the system detects a load change and reinitializes the voltage reference to _V<sub>s0</sub>_ .  
5. The emulation process restart, allowing the system to operate again as a PV emulator.






