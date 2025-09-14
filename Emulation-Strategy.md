The PV emulator proposed in this repository relies heavily on the understanding of the **current–voltage (I–V) characteristic curve** of photovoltaic (PV) modules. For this reason, it is important to examine this curve in more detail.

## _Current–Voltage (I–V) Characteristic Curve_  

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

- **Simplified exponential model**  
- **Single-diode model**

Regardless of which of these two models is considered, both are capable of reproducing the complete I–V characteristic curve of a specific PV module using only a few parameters that are easily available in the datasheet.

## _Emulation Strategy_  



<img width="500" height="400" alt="Emulator_parameters" src="https://github.com/user-attachments/assets/55ae40e4-b984-4d14-a78b-da1e50bd402b" />

