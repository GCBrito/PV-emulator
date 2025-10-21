import matplotlib.pyplot as plt
import numpy as np

# Configura a fonte padrão para 'Times New Roman', similar ao comando do MATLAB
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# Equivalente a 'clear, clc, close all' do MATLAB
# (limpar variáveis não é necessário em um script, 'clc' não tem equivalente direto)
plt.close('all')

# %% 1) Display of the three curves + experimental points
# Cria a figura com cor de fundo branca e tamanho específico (convertido de pixels para polegadas)
fig = plt.figure(facecolor='w', figsize=(9, 5), dpi=100)
ax = fig.add_subplot(1, 1, 1)

# Ativa o grid, equivalente a 'hold on; grid on;' (hold on é o padrão no Matplotlib)
ax.grid(True)

# Single diode emulator - KC85TS CS6P-250P (S = 765 W/m^2 and T = 44.5°C)
# V_emulator = np.array([
#     0,
#     22.80755,
#     24.17536,
#     25.77529,
#     27.24029,
#     28.71549,
#     29.93957,
#     30.82569,
#     32.12266,
#     33.61550,
#     33.92,
# ])
# 
# I_emulator = np.array([
#     6.79525,
#     6.728934,
#     6.696453,
#     6.513734,
#     6.299920,
#     5.674942,
#     5.156351,
#     4.187371,
#     2.577053,
#     0.723580,
#     0,
# ])

# Single diode emulator - KC85TS CS6P-250P (S = 556 W/m^2 and T = 33°C)
V_emulator = np.array([
    0,
    15.99633,
    21.38497,
    24.11979,
    25.47998,
    27.13490,
    28.61119,
    30.10375,
    31.92938,
    33.62070,
    34.34491,
    34.9,
])
I_emulator = np.array([
    4.93461,
    4.926764,
    4.906911,
    4.890384,
    4.836957,
    4.766093,
    4.515424,
    4.260068,
    3.011903,
    1.829092,
    0.796171,
    0,
])

# plot single diode emulator
h1, = ax.plot(V_emulator, I_emulator, '-b', linewidth=2)

# Real panel data - KC85TS CS6P-250P (S = 765 W/m^2 and T = 44.5°C)
# V_pannel = np.array([
#     0.047996217771366645,
#     2.5473136313978317,
#     4.686489351002983,
#     7.605990366562553,
#     10.181947580363017,
#     12.610309040298777,
#     14.79342630980413,
#     17.20541669131795,
#     19.415489834424953,
#     21.816632013769198,
#     24.03214755291164,
#     26.237334370932288,
#     28.26322957001055,
#     29.749962170030145,
#     31.046437882063465,
#     31.92801141798128,
#     32.689215781035756,
#     33.25368190774834,
#     33.59447586154349,
#     34.38443600476654
# ])
# 
# I_pannel = np.array([
#     7.007823910286788,
#     7.0079635369891475,
#     7.001379133258391,
#     7.003776871281191,
#     6.93017776227956,
#     6.9347826983712295,
#     6.841049913320265,
#     6.845653933494394,
#     6.8524813120762555,
#     6.83697299652661,
#     6.848269952562965,
#     6.68079538185797,
#     6.2719700120743145,
#     5.5323882664885256,
#     4.567097566839892,
#     3.6531803358379182,
#     2.83534576727134,
#     2.1046510438953376,
#     1.3538320938536454,
#     0.08683712275429656
# ])

# Real panel data - KC85TS CS6P-250P (S = 556 W/m^2 and T = 33°C)
V_pannel = np.array([
   0.11870949058476832,
   2.503084084256967,
   5.175401505238524,
   7.590304253402278,
   9.815231816381960,
   12.187280070347374,
   14.589991865444773,
   16.821065673119850,
   19.028722204632890,
   21.272139279014567,
   23.472448005559684,
   25.672959815744026,
   28.063920664439690,
   29.917176441633483,
   31.601112724298030,
   32.722065096396480,
   33.616384803786715,
   34.327349424005450,
   34.741398772039830,
   34.904501246242180,
])
I_pannel = np.array([
    5.110461288554346,
    5.100720374022223,
    5.111095915194137,
    5.118929111074918,
    5.121718383239790,
    5.122016080125351,
    5.119807586122182,
    5.120087590427699,
    4.957212044814572,
    4.944943397608992,
    4.962789821397647,
    4.950515789494053,
    4.873004614033777,
    4.546931980897691,
    3.879472637469006,
    3.076400466372332,
    2.250709493400178,
    1.347184263076274,
    0.845228194983978,
    0.290529789361583,
])

# plot real panel data
h2, = ax.plot(V_pannel, I_pannel, 'ko', 
             markersize=6, markerfacecolor='k')

ax.set_xlabel('Voltage (V)', fontsize=30)
ax.set_ylabel('Current (A)', fontsize=30)

# Axis limits to match the provided figure
ax.set_xlim(0, 36)
ax.set_ylim(0, 5.5)

# Legenda
ax.legend([h1, h2], 
          ['PV Emulator – Single Diode Model', 'Experimental Results – Real PV Panel'], 
          loc='best', fontsize=20) # Ajustei o fontsize da legenda para melhor visualização

# Axis ticks font configuration
ax.tick_params(axis='both', which='major', labelsize=30)

# Mostra o gráfico
plt.tight_layout() # Ajusta o padding para que os labels não fiquem cortados
plt.show()