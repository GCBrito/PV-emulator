import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

plt.close('all')

## 1) Display of the three curves + experimental points
# figure('Color','w','Position',[100 100 900 500]);
fig, ax = plt.subplots(figsize=(12, 7), facecolor='w')
ax.grid(True)

# Single diode emulator - CanadianSolar CS6P-250P (S = 765 W/m^2 and T = 44.5°C)
# V_emulator = np.array([0, 22.80755, 24.17536, 25.77529, 27.24029, 28.71549, 29.93957, 30.82569, 32.12266, 33.61550, 33.92])
# I_emulator = np.array([6.79525, 6.728934, 6.696453, 6.513734, 6.299920, 5.674942, 5.156351, 4.187371, 2.577053, 0.723580, 0])

# Single diode emulator - CanadianSolar CS6P-250P (S = 556 W/m^2 and T = 33°C)
# V_emulator = np.array([0, 15.99633, 21.38497, 24.11979, 25.47998, 27.13490, 28.61119, 30.10375, 31.92938, 33.62070, 34.34491, 34.9])
# I_emulator = np.array([4.93461, 4.926764, 4.906911, 4.890384, 4.836957, 4.766093, 4.515424, 4.260068, 3.011903, 1.829092, 0.796171, 0])

# Single diode emulator - Vertex N TSM-NEG21C.20 695W (S = 200 W/m^2 and T = 25°C)
# V_emulator = np.array([0, 17.31385, 29.14563, 33.30596, 36.31013, 38.27775, 40.10532, 41.2025, 42.82103, 44.22889, 44.3])
# I_emulator = np.array([3.66, 3.650519, 3.634901, 3.602349, 3.529229, 3.277086, 3.042891, 2.424462, 1.360115, 0.434306, 0])

# Single diode emulator - Kyocera KC200GT (S = 511 W.m^-2 and T = 54.3 ºC)
V_emulator = np.array([0, 11.90924, 15.66191, 18.63553, 21.14057, 22.63848, 23.60215, 24.54963, 25.56798, 26.31991, 26.9799, 27.3])
I_emulator = np.array([4.2, 4.188322, 4.172224, 4.050412, 3.926395, 3.60058, 3.390971, 2.825413, 2.187756, 1.699926, 0.723199, 0])

h1, = ax.plot(V_emulator, I_emulator, '-b', linewidth=2, label='PV Emulator – Single Diode Model')

# Real panel data - CanadianSolar CS6P-250P (S = 765 W/m^2 and T = 44.5°C)
# V_pannel = np.array([0.047996, 2.54731, 4.68648, 7.60599, 10.1819, 12.6103, 14.7934, 17.2054, 19.4154, 21.8166, 24.0321, 26.2373, 28.2632, 29.7499, 31.0464, 31.9280, 32.6892, 33.2536, 33.5944, 34.3844])
# I_pannel = np.array([7.0078, 7.0079, 7.0013, 7.0037, 6.9301, 6.9347, 6.8410, 6.8456, 6.8524, 6.8369, 6.8482, 6.6807, 6.2719, 5.5323, 4.5670, 3.6531, 2.8353, 2.1046, 1.3538, 0.0868])

# Real panel data - CanadianSolar CS6P-250P (S = 556 W/m^2 and T = 33°C)
# V_pannel = np.array([0.1187, 2.5030, 5.1754, 7.5903, 9.8152, 12.1872, 14.5899, 16.8210, 19.0287, 21.2721, 23.4724, 25.6729, 28.0639, 29.9171, 31.6011, 32.7220, 33.6163, 34.3273, 34.7413, 34.9045])
# I_pannel = np.array([5.1104, 5.1007, 5.1110, 5.1189, 5.1217, 5.1220, 5.1198, 5.1200, 4.9572, 4.9449, 4.9627, 4.9505, 4.8730, 4.5469, 3.8794, 3.0764, 2.2507, 1.3471, 0.8452, 0.2905])

# Real panel data - Vertex N TSM-NEG21C.20 695W (S = 200 W/m^2 and T = 25°C)
# V_pannel = np.array([0.1678, 11.8815, 21.8413, 30.5910, 32.8237, 34.4093, 35.8007, 37.2501, 38.8028, 40.5299, 41.9849, 43.2068, 44.4279])
# I_pannel = np.array([3.7614, 3.7416, 3.7507, 3.7431, 3.7516, 3.7559, 3.7186, 3.5650, 3.2617, 2.7423, 2.0733, 1.2878, 0.0743])

# Single diode emulator - Kyocera KC200GT (S = 511 W.m^-2 and T = 54.3 ºC)
V_pannel = np.array([0.0663, 2.7458, 4.7501, 6.7383, 8.5551, 10.3986, 12.2473, 14.0855, 15.9289, 17.7671, 19.4713, 21.3199, 22.8683, 24.2877, 25.4336, 26.2901, 26.8304, 27.4188, 27.8465, 28.2476])
I_pannel = np.array([4.1435, 4.1571, 4.1616, 4.1618, 4.1552, 4.1620, 4.0997, 4.0909, 4.0331, 4.0466, 4.0445, 3.9221, 3.6816, 3.2696, 2.7439, 2.2115, 1.6769, 1.0888, 0.5074, 0.0574])

h2 = ax.scatter(V_pannel, I_pannel, c='k', marker='o', s=60, edgecolors='k', label='Experimental Results – Real PV Panel')

ax.set_xlabel('Voltage (V)', fontsize=30)
ax.set_ylabel('Current (A)', fontsize=30)

ax.legend(loc='best', fontsize=12, frameon=True)

ax.tick_params(axis='both', which='major', labelsize=30)

for item in ([ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontname('Times New Roman')

plt.tight_layout()
plt.show()