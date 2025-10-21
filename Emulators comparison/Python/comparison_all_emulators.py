import numpy as np
import matplotlib.pyplot as plt

# Script Python : comparaison des 3 courbes + points expérimentaux

# Équivalent de 'clear; clc; close all' en MATLAB
plt.close('all')

# Appliquer la police 'Times New Roman' à tous les éléments de texte
# C'est l'équivalent des commandes 'set(findall...)' de MATLAB
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# %% 1) Affichage des trois courbes + points expérimentaux
# Créer la figure avec une couleur de fond blanche et une taille spécifique
fig, ax = plt.subplots(figsize=(9, 5), dpi=100, facecolor='w')

# Activer la grille et permettre de superposer les graphiques ('hold on' est le défaut)
ax.grid(True)

# % Émulateur complet
# V_complet = np.array([
#     0,
#     13.13,
#     18.31,
#     20.56,
#     21.88,
#     22.68,
#     23.52,
#     24.61,
#     25.46,
#     26.1,
#     26.69,
#     27.24,
#     27.52,
#     27.75,
#     28.18,
#    ])
#
# I_complet = np.array([
#     4.12,
#     4.03,
#     3.96,
#     3.87,
#     3.72,
#     3.63,
#     3.38,
#     3,
#     2.7,
#     2.19,
#     1.72,
#     0.99,
#     0.59,
#     0.26,
#     0,
# ])

# % Émulateur complet - KC85TS
# V_complete = np.array([
#     0,
#     13.20,
#     14.49,
#     15.37,
#     16.52,
#     17.29,
#     18.48,
#     19.51,
#     20.11,
#     20.73,
#     21.23,
#     21.55,
#     21.7
# ])
#
# I_complete = np.array([
#     5.34,
#     5.33,
#     5.31,
#     5.28,
#     5.15,
#     5.04,
#     4.29,
#     3.60,
#     2.77,
#     1.85,
#     0.89,
#     0.29,
#     0
# ])

# Émulateur complet - Uni‑Solar ES‑62T (S = 1000 W/m^2 et T = 25°C)
V_complete = np.array([
    0,
    11.90,
    12.39,
    13.38,
    14.47,
    15.04,
    16.60,
    18.09,
    18.88,
    19.93,
    20.58,
    21,
])
I_complete = np.array([
    5.1,
    4.62,
    4.56,
    4.44,
    4.21,
    4.08,
    3.38,
    2.68,
    2.09,
    1.15,
    0.46,
    0,
])

h1, = ax.plot(V_complete,  I_complete,   '-b',  linewidth=2)

# % Émulateur simple 
# V_simple = np.array([
#     0,
#     15.59,
#     17.23,
#     18.55,
#     19.60,
#     20.42,
#     21.61,
#     22.88,
#     24.18,
#     25.40,
#     26.10,
#     26.69,
#     27.49,
#     27.97,
#     28.18,
# ])
#
# I_simple = np.array([
#     4.12,
#     4.08,
#     4.06,
#     4.03,
#     3.98,
#     3.94,
#     3.80,
#     3.63,
#     3.18,
#     2.76,
#     2.27,
#     1.84,
#     0.87,
#     0.26,
#     0,
# ])

# % Émulateur simple - KC85TS
# V_simple = np.array([
#     0,
#     13.24,
#     14.25,
#     15.19,
#     15.95,
#     17.23,
#     18.29,
#     19.49,
#     20.09,
#     20.62,
#     21.11,
#     21.58,
#     21.7
# ])
#
# I_simple = np.array([
#     5.34,
#     5.31,
#     5.29,
#     5.26,
#     5.20,
#     5.04,
#     4.61,
#     4.06,
#     3.36,
#     2.71,
#     1.48,
#     0.29,
#     0
# ])

# Émulateur simple -Uni‑Solar ES‑62T (S = 1000 W/m^2 et T = 25°C)
V_simple = np.array([
    0,
    10.14,
    11.60,
    12.88,
    14.14,
    14.96,
    16.69,
    17.92,
    18.90,
    19.91,
    20.49,
    21,
])
I_simple = np.array([
    5.1,
    4.81,
    4.69,
    4.53,
    4.29,
    4.11,
    3.39,
    2.87,
    2.22,
    1.31,
    0.61,
    0,
])

h2, = ax.plot(V_simple, I_simple, '-r', linewidth=2)  # plot modèle exponentielle simple

# % Émulateur commercial
# V_commercial = np.array([
#     0,
#     1.43,
#     4.82,
#     9.5,
#     15.9,
#     18.72,
#     19.9,
#     20.51,
#     22.06,
#     23.07,
#     23.95,
#     24.58,
#     24.79,
#     25.39,
#     25.90,
#     26.43,
#     27.00,
#     27.57,
#     28.09,
#     28.18,
#    ])
#
# I_commercial = np.array([
#     4.12,
#     4.11,
#     4.09,
#     4.05,
#     3.97,
#     3.89,
#     3.85,
#     3.83,
#     3.73,
#     3.61,
#     3.37,
#     2.96,
#     2.79,
#     2.3,
#     1.87,
#     1.44,
#     0.97,
#     0.5,
#     0.07,
#     0,
# ])

# Émulateur commercial - Uni‑Solar ES‑62T (S = 1000 W/m^2 et T = 25°C)
V_commercial = np.array([
    0, 3.05, 5.95, 7.42, 10.72, 11.25, 12.2, 12.93, 13.54, 14.22,
    14.28, 14.6, 14.9, 15.08, 15.31, 15.63, 15.84, 16.05, 16.26, 16.46,
    16.71, 17.02, 17.31, 17.54, 17.7, 18, 18.54, 19.01, 19.56, 20.02,
    20.5, 20.94, 21,
])
I_commercial = np.array([
    5.1, 5.00, 4.90, 4.83, 4.62, 4.64, 4.49, 4.40, 4.32, 4.22,
    4.21, 4.18, 4.10, 4.06, 4.02, 3.95, 3.88, 3.81, 3.75, 3.69,
    3.58, 3.44, 3.29, 3.16, 3.049, 2.85, 2.45, 2.08, 1.63, 1.26,
    0.85, 0.05, 0,
])

h3, = ax.plot(V_commercial, I_commercial, '-m', linewidth=2)  # plot émulateur commercial

# % Données panneau réel
# tension_orig = np.array([...])
# courant_orig = np.array([...])

# % Données panneau réel - KC85TS
# V_pannel = np.array([...])
# I_pannel = np.array([...])

# Données panneau réel - Uni‑Solar ES‑62T (S = 1000 W/m^2 et T = 25°C)
V_pannel = np.array([
    0.02882, 0.89345, 1.90218, 2.78603, 3.88122, 5.03406, 6.05240, 
    7.04192, 7.99301, 9.14585, 10.23144, 11.21135, 11.84541, 12.54672, 
    13.10393, 13.58428, 14.08384, 14.53537, 14.92926, 15.26550, 
    15.66900, 16.03406, 16.37991, 16.76419, 17.17729, 17.61921, 
    18.13799, 18.65677, 19.09869, 19.63668, 20.16507, 20.57817, 20.94323
])
I_pannel = np.array([
    5.06297, 5.02784, 4.98919, 4.95405, 4.90838, 4.86270, 4.82054, 
    4.78189, 4.72919, 4.67297, 4.61676, 4.56757, 4.53243, 4.49730, 
    4.44811, 4.39189, 4.32865, 4.24784, 4.15297, 4.05811, 3.92108, 
    3.77351, 3.59784, 3.39405, 3.12703, 2.84595, 2.46297, 2.04838, 
    1.67243, 1.19108, 0.73784, 0.34784, 0.00703
])

# plot données panneau réel
h4, = ax.plot(V_pannel, I_pannel, 'ko', markersize=6, markerfacecolor='k')

ax.set_xlabel('Tension (V)', fontsize=30)
ax.set_ylabel('Courant (A)', fontsize=30)
ax.legend(handles=[h1, h2, h3, h4],
          labels=['Émulateur open-source modèle à une diode', 'Émulateur open-source modèle exp. simplifié', 'Émulateur commercial', 'Données panneau réel'],
          loc='best', fontsize=14) # Ajustement de la taille de la police pour la lisibilité

# Configurations de police pour les graduations des axes
ax.tick_params(axis='both', which='major', labelsize=30)

# plt.tight_layout() # Dé-commentez si les labels sont coupés
plt.show()