#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 18:49:49 2023

@author: kevin
"""
import numpy as np
import matplotlib.pyplot as plt

# initialising the resistors and inductors
R = 2  # series resistance
R1 = 1  # resistance 1
R2 = 2  # resistance 2
L1 = 1e-3  # inductance 1
L2 = 0.5e-3  # inductance 2
IL0 = 0  # initial condition 1
IL1 = 0  # initial condition 2

# Defining the voltage source
T = 6e-3  # simulation time
dt = 1e-8  # time step
V0 = 1  # voltage amplitude for source

t = np.arange(0, T, dt)  # creating the time array
K = int(len(t)/2)
# cosine voltage
# V=np.cos(2*np.pi*1e3*t)
V = np.concatenate((np.ones(K), -np.ones(K)))
I1 = np.zeros(len(t)+1, dtype=np.float32)
I2 = np.zeros(len(t)+1, dtype=np.float32)
Vinst = np.zeros(len(t), dtype=np.float32)
I1[0], I2[0] = IL0, IL1
# main simulation
for i in range(len(t)):
    Vinst[i] = V[i]-(I1[i]+I2[i])*R
    dI1, dI2 = (Vinst[i]-R1*I1[i])/(R1+L1/dt), (Vinst[i]-R2*I2[i])/(R2+L2/dt)
    I1[i+1] += dI1+I1[i]
    I2[i+1] += dI2+I2[i]
plt.subplot(4, 1, 1)
plt.plot(t, I1[1:], color='black')
plt.subplot(4, 1, 2)
plt.plot(t, Vinst)
plt.subplot(4, 1, 3)
plt.plot(t, I2[1:])
plt.subplot(4, 1, 4)
plt.plot(t, V)
plt.show()
