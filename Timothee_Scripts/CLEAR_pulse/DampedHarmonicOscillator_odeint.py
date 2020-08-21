# -*- coding: utf-8 -*-
"""
@author: Timothee Guerra

Solve the differential equation of a damped harmonic oscillator using odeint

"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint
    

def CLEARDrivenForce(t, ti, tf, P, relAmp, Tsegm_up, Tseg_down):
    if t < ti:
        return 0
    elif t < ti + Tsegm_up/2:
        return np.sqrt(P)*relAmp[0]
    elif t < ti + Tsegm_up:
        return np.sqrt(P)*relAmp[1]
    elif t < tf:
        return np.sqrt(P)
    elif t < tf + Tseg_down/2:
        return -np.sqrt(P)*relAmp[2]
    elif t < tf + Tseg_down:
        return np.sqrt(P)*relAmp[3]
    else:
        return 0
    

def squareDrivenForce(t, ti, tf, P):
    if t<ti:
        return 0
    elif t<tf:
        return np.sqrt(P)
    else:
        return 0


def dynamics(Y, t, G, Wo, func, funcArgs):
        
    x, v = Y[0], Y[1]
    
    dx_dt = v
    dv_dt = -2*G*v - Wo**2*x + func(t, *funcArgs)
    
    return np.array( [dx_dt, dv_dt] )


IC = [0, 0] # the initial condition, [ x, v] at t=0
ti, Dt = 10, 100
tf = ti+Dt
Tsegm_up, Tseg_down = 20, 20
relAmp = (1.1, 0.9, 0.2, 0.2)

P = 1
G, Wo = 0.1, np.sqrt(2)

T = np.linspace(ti-10, ti+Dt+50, 5000)

# paramfuncSquare = [ti, tf, P]
# simuParamSquare = (G, Wo, squareDrivenForce, paramfuncSquare)
# YmatSquare = odeint(dynamics, IC, T, args=simuParamSquare)

paramfuncClear = [ti, tf, P, relAmp, Tsegm_up, Tseg_down]
simuParamCLEAR = (G, Wo, CLEARDrivenForce, paramfuncClear)
YmatCLEAR = odeint(dynamics, IC, T, args=simuParamCLEAR)

steadyStateUp = np.sqrt(P)/Wo**2
steadyStateDown = 0

# listResult = [YmatSquare, YmatCLEAR]
# listPrint = ['Square', 'CLEAR']
# for result, dprint in zip(listResult, listPrint):
#     popTime, depopTime = None, None
#     for i in range(len(result[:, 0])):
#         indTf = np.argwhere(T==tf)[0,0]
#         if np.allclose(result[i:indTf, 0], steadyStateUp, rtol=1e-3, atol=0) and popTime is None:
#             popTime = T[i]
#         if np.allclose(result[i:None, 0], steadyStateDown, atol=1e-3) and depopTime is None:
#             depopTime = T[i]

#     print('{} pulse : populating time {}ns - depopulating time {}ns'.format(dprint, popTime, depopTime-tf))



#####################################
#####----  Plot solutions ------#####
#####################################


fig = plt.figure(num=1, figsize=[18,8])
ax1 = fig.add_subplot(1, 1, 1)

color = 'tab:blue'

# ax1.plot(T, YmatSquare[:, 0], ls=':', color=color, label='Square pulse')
ax1.plot(T, YmatCLEAR[:, 0], ls='-', color=color, label='CLEAR pulse')
ax1.set_xlabel('time (ns)')
ax1.set_ylabel('Cavity population', color=color)
ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'grey'
# ax1.plot(T, [squareDrivenForce(t, *paramfuncSquare) for t in T], ls='--', color=color)
ax1.plot(T, [CLEARDrivenForce(t, *paramfuncClear) for t in T], ls='-', color=color)

# ax2.set_ylabel('Readout shape', color=color)  # we already handled the x-label with ax1
# ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()