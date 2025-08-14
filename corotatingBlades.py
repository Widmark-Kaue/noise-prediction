#%% imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from src.noise import farfield

dir = Path().joinpath('validate_data', 'corotatingBlades')


#%% load mics
# Array of mics
x, z, y = np.loadtxt(dir.joinpath('array24mics.txt'), unpack=True)
y = np.full_like(y, 1.06)

mic = np.array([x, y, z]).T
S = np.sqrt(x**2 + y**2 + z**2)

# Propeller
D = 12 * 0.0254
B = 2
phi = np.linspace(0, 2*np.pi, 25)
reff = D/2 *0.8 
r = D/2

y1 = reff*np.cos(phi)
z1 = reff*np.sin(phi)
x1 = np.zeros(len(y1))

#%% Plot
# Reference lines
n = 100
line_x = np.array([
    np.linspace(0, 1.06, n),
    np.zeros(n),
    np.zeros(n)
])

line_y = np.array([
    np.zeros(n),
    np.linspace(0, 1.06, n),
    np.zeros(n)
])
line_z = np.array([
    np.zeros(n),
    np.zeros(n),
    np.linspace(0, 1.06, n),
])
ax = plt.figure().add_subplot(projection='3d')
ax.plot(x1, y1, z1, color='b')
ax.plot(x, y, z, 'ko')
ax.plot(line_x[0, :], line_x[1, :], line_x[2, :], 'k--')
ax.plot(line_y[0, :], line_y[1, :], line_y[2, :], 'k--')
ax.plot(line_z[0, :], line_z[1, :], line_z[2, :], 'k--')

ax.set_xlim(-0.5, 1.06)
ax.set_zlim(-0.5, 1.06)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# ax.view_init(elev=15., azim=25, roll=0)
plt.tight_layout()
plt.show(block = False)

#%% Plot array mics 2D
plt.figure()
plt.plot(x, z, 'ko')
plt.axhline(0, color = 'k', linestyle = '--')
plt.axvline(0, color = 'k', linestyle = '--')
plt.axis('equal')

plt.xlabel('x')
plt.ylabel('z')
plt.show(block = False)

#%% Noise computation
Pref = 2e-5         # Pa

data = pd.read_excel(dir.joinpath('corotating_data.xlsx'))
case = farfield(microphones=mic)

T = data['Thrust [N]']
Q = data['Torque [N.m]']
Mrot = data['RPM'] * 2 * np.pi * r/ 60 / case.sound_speed
Mx = 0

T_apc = 5.501
Q_apc = 0.112

prms = case.garrickWatkinsReff(
        number_of_harmonics=1,
        number_of_blades=B, 
        reff=reff,
        # loading=[T[0], Q[0]],
        # Mrot=Mrot[0],
        loading=[T_apc, Q_apc],
        Mrot = 0.19 * 0.8, 
        Mx=Mx
    )

SPL = 20*np.log10(prms/Pref)
print(SPL)

ax = plt.figure().add_subplot(projection='3d')
ax.plot(x, z, SPL.T,  'ko')
plt.show()


# %% Nois
