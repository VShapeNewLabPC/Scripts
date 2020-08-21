# -*- coding: utf-8 -*-
"""
@author: Timothee Guerra

Solve the differential equation of a damped harmonic oscillator using odeint

"""



import time
import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize
from scipy.integrate import odeint


MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



def instantRespStepUp(t, ti, P, phiParams):
    w = phiParams['omega']
    B = phiParams['beta']
    y0 = phiParams['y(0)']
    dy0 = phiParams['dy(0)']
    
    if t < ti:
        return y0, dy0
    
    if P == 0:
        return y0, dy0
    
    O = np.sqrt(w**2 - B**2)
    tp = t - ti
    
    A1 = -y0*O*w**2/P
    A2 = B*A1/O - dy0*w**2/P
        
    osci = (O + A1)*np.cos(O*tp) + (B + A2)*np.sin(O*tp)
    cte = P/w**2
    yt = cte*(1 - np.exp(-B*tp)*osci/O) 
    
    dosci = -(O+A1)*np.sin(O*tp) + (B+A2)*np.cos(O*tp)
    dyt = (cte*B*osci/O - cte*dosci)*np.exp(-B*tp)
    
    return yt, dyt

def instantRespStepDown(t, ti, P, phiParams):
    w = phiParams['omega']
    B = phiParams['beta']
    y0 = phiParams['y(0)']
    dy0 = phiParams['dy(0)']
    
    if t < ti:
        return y0, dy0
    
    if P == 0:
        return y0, dy0
    
    O = np.sqrt(w**2 - B**2)
    tp = t - ti
    
    A1 = y0*O*w**2/P - O
    A2 = B*A1/O + dy0*w**2/P
        
    osci = (O + A1)*np.cos(O*tp) + (B + A2)*np.sin(O*tp)
    cte = P/w**2
    yt = cte*np.exp(-B*tp)*osci/O
    
    dosci = -(O+A1)*np.sin(O*tp) + (B+A2)*np.cos(O*tp)
    dyt = (-cte*B*osci/O + cte*dosci)*np.exp(-B*tp)
    
    return yt, dyt



class pulseShaping(object):
    
    phiParams = None
    ratedPower = None
    
    timeStart = None
    timeStop = None
    timeVec = None

    pulseStart = None
    pulseEnd = None
    pulseDuration = None
    
    upDuration = None
    upPower = None
    
    downDuration = None
    downPower = None
    
    drivingPulse = None
    responsePulse = None
    
    
    def __init__(self, phiParams, pulseParam, ti, tf, upParam=None, downParam=None):
        self.phiParams = phiParams

        self.pulseStart = pulseParam['start']
        self.pulseDuration = pulseParam['duration']
        self.pulseEnd = self.pulseStart + self.pulseDuration
        
        self.ratedPower = pulseParam['power']
        
        if upParam is not None:
            if type(upParam) == dict:
                self.upDuration = upParam['duration']
                self.upPower = upParam['power']
            else:
                nbUpParam = len(upParam)
                self.upDuration = list(upParam[:nbUpParam//2])
                self.upPower = list(upParam[nbUpParam//2:])
            
        if downParam is not None:
            if type(downParam) == dict:
                self.downDuration = downParam['duration']
                self.downPower = downParam['power']
            else:
                nbDownParam = len(downParam)
                self.downDuration = list(downParam[:nbDownParam//2])
                self.downPower = list(downParam[nbDownParam//2:])
            
                
        self.timeStart = ti
        self.timeStop = tf
        self.timeVec = np.linspace(ti, tf, 5000)
        
    def drivingStepUpPulse(self):
        self.drivingPulse = [self.ratedPower if t>self.pulseStart else 0 
                             for t in self.timeVec]

    def drivingStepDownPulse(self):
        self.drivingPulse = [self.ratedPower if t<self.pulseStart else 0 
                             for t in self.timeVec]

    def drivingSquarePulse(self):
        stepUpStart = np.array([self.ratedPower if t>self.pulseStart else 0 
                                for t in self.timeVec])
        stepUpEnd = np.array([self.ratedPower if t>self.pulseStart+self.pulseDuration else 0 
                                for t in self.timeVec])
        
        self.drivingPulse = stepUpStart - stepUpEnd
        
    def drivingClearPartPulse(self, typePart):        
        
        power, duration = self.getPowerDuration_forClear(typePart)
        self.drivingPulse = np.zeros(np.shape(self.timeVec), dtype=float)
        for i, t in enumerate(self.timeVec):
            
            limitTime = 0
            for RP, dt in zip(power, duration):
                
                limitTime += dt
                if t<=limitTime:
                    self.drivingPulse[i] = RP * self.ratedPower
                    break
    
    
    def responseStepUp(self):
        ti, P, phiParams = self.pulseStart, self.ratedPower, self.phiParams
        self.responsePulse = [instantRespStepUp(t, ti, P, phiParams)[0] for t in self.timeVec]
    
    def responseStepDown(self):
        ti, P, phiParams = self.pulseStart, self.ratedPower, self.phiParams
        self.responsePulse = [instantRespStepDown(t, ti, P, phiParams)[0] for t in self.timeVec]

    def responseSquare(self):
        y0, dy0 = None, None
        
        self.responsePulse = np.full(np.shape(self.timeVec), self.phiParams['y(0)'], dtype=float)
        for i, t in enumerate(self.timeVec):
            if t < self.pulseStart:
                continue
            
            if t < self.pulseEnd:
                ti, P, phiParams = self.pulseStart, self.ratedPower, self.phiParams.copy()
                y0, dy0 = instantRespStepUp(t, ti, P, phiParams)
                self.responsePulse[i] = y0
            else:
                ti, P, phiParams = self.pulseEnd, self.ratedPower, self.phiParams.copy()
                phiParams['y(0)'], phiParams['dy(0)'] = y0, dy0
                self.responsePulse[i] = instantRespStepDown(t, ti, P, phiParams)[0]
     
    def getPowerDuration_forClear(self, typeClear):
        if typeClear == 'up':
            power = [0.0] + self.upPower + [1.0]
            duration = [self.pulseStart] + self.upDuration
            duration = duration + [self.timeStop-np.sum(duration)]      
        elif typeClear == 'down':
            power = [1.0] + self.downPower + [0.0]
            duration = [self.pulseStart] + self.downDuration
            duration = duration + [self.timeStop-np.sum(duration)]
        elif typeClear == 'all':
            power = [0.0] + self.upPower + [1.0] + self.downPower + [0.0]
            duration = [self.pulseStart]
            duration += self.upDuration 
            duration += [self.pulseEnd-np.sum(duration)]
            duration += self.downDuration
            duration += [self.timeStop-np.sum(duration)]
        else:
            power, duration = None, None
            
        return power, duration   
     
    
    def responseClearPart(self, typePart):   
        power, duration = self.getPowerDuration_forClear(typePart)
        nbSegment = len(duration)
        
        IC = [{'y(end)': self.phiParams['y(0)'], 
               'dy(end)': self.phiParams['dy(0)']} for _ in range(nbSegment)]
        
        self.responsePulse = np.full(np.shape(self.timeVec), self.phiParams['y(0)'], dtype=float)
        for i, t in enumerate(self.timeVec):
            limitTime = self.pulseStart
            if t < limitTime:
                continue
            
            for numSegment in range(1, nbSegment):
                startTime = limitTime
                limitTime += duration[numSegment]
                if t <= limitTime: 
                    ti, P, phiParams = startTime, power[numSegment]*self.ratedPower, self.phiParams.copy()
                    phiParams['y(0)'], phiParams['dy(0)'] = IC[numSegment-1]['y(end)'], IC[numSegment-1]['dy(end)']
                    
                    if numSegment == nbSegment-1 and typePart != 'up':
                        P = power[numSegment-1]*self.ratedPower
                        y0, dy0 = instantRespStepDown(t, ti, P, phiParams)
                        IC[numSegment]['y(end)'], IC[numSegment]['dy(end)'] = y0, dy0
                        self.responsePulse[i] = y0
                    else:
                        y0, dy0 = instantRespStepUp(t, ti, P, phiParams)
                        IC[numSegment]['y(end)'], IC[numSegment]['dy(end)'] = y0, dy0
                        self.responsePulse[i] = y0
    
    
    def getPopulationTime(self, tol=1e-3):
        steadyStateUp = self.ratedPower/self.phiParams['omega']**2
        
        T, ti, tf = self.timeVec, self.pulseStart, self.pulseEnd
        indTi = np.argwhere(T==T[T>=ti][0])[0,0]
        indTf = np.argwhere(T==T[T>=tf][0])[0,0]
        for i in range(indTi, len(T)):
            if np.allclose(self.responsePulse[i:indTf], steadyStateUp, rtol=tol, atol=0):
                return T[i]
        
        return 1e15
    
    def getDepopulationTime(self, tol=1e-3):
        steadyStateDown = 0
        
        T, ti = self.timeVec, self.pulseStart
        indTi = np.argwhere(T==T[T>=ti][0])[0,0]
        for i in range(indTi, len(T)):
            if np.allclose(self.responsePulse[i:], steadyStateDown, atol=tol):
                return T[i]
        
        return 1e15
    
    
    
     
def minimumPopulationTime(upParam, phiParam, pulseParam, ti, tf, tol):

    pulse = pulseShaping(phiParam, pulseParam, ti, tf, upParam=upParam)
    pulse.responseClearPart('up')
    popTime = pulse.getPopulationTime(tol)
    print('\rPopulation duration found :   {:.3f}'.format(popTime), end='')

    return popTime
    

def getOptimizeClearUpParam(upParam, phiParam, pulseParam, ti, tf, tol=1e-3):
    duration = upParam['duration']
    power = upParam['power']
    x0 = [*duration, *power]
    args = (phiParam, pulseParam, ti, tf, tol)
    
    startUp = time.time()
    resultUp = optimize.minimize(minimumPopulationTime, x0, args=args, method='Nelder-Mead')
    endUp = time.time()
    print('\nDuration minimization of clear up : {}s'.format(endUp-startUp))
    return resultUp.x


def minimumDepopulationTime(downParam, phiParam, pulseParam, ti, tf, tol):
    backup = phiParam['y(0)']
    phiParam['y(0)'] = pulseParam['power'] / phiParam['omega']**2
    pulse = pulseShaping(phiParam, pulseParam, ti, tf, downParam=downParam)
    pulse.responseClearPart('down')
    popTime = pulse.getDepopulationTime(tol)
    print('\rDepopulation duration found :   {:.3f}'.format(popTime), end='')
    
    phiParam['y(0)'] = backup

    return popTime
   

def getOptimizeClearDownParam(upParam, phiParam, pulseParam, ti, tf, tol=1e-3):
    duration = downParam['duration']
    power = downParam['power']
    x0 = [*duration, *power]
    args = (phiParam, pulseParam, ti, tf, tol)
    
    startDown = time.time()
    resultDown = optimize.minimize(minimumDepopulationTime, x0, args=args, method='Nelder-Mead')
    endDown = time.time()
    print('\nDuration minimization of clear down : {}s'.format(endDown-startDown))
    return resultDown.x
     
##################################################
############------- PARAMETERS -------############
##################################################


ti, tf = 0, 200

phiParam = {'omega': np.sqrt(2), 'beta': 0.2, 'y(0)': 0, 'dy(0)': 0}
pulseParam = {'start': 10, 'duration': 100, 'power': 1}
upParam = {'duration': [10, 10], 'power': [1.1, 0.9]}
downParam = {'duration': [10, 10], 'power': [-0.2, 0.2]}

tol = 1e-3

############################################
############------- MAIN -------############
############################################


popTimeClear, depopTimeClear = None, None
popTimeSquare, depopTimeSquare = None, None
pulseSquare, pulseClear = None, None

# upParam = getOptimizeClearUpParam(upParam, phiParam, pulseParam, ti, tf, tol=tol)
# downParam = getOptimizeClearDownParam(upParam, phiParam, pulseParam, ti, tf, tol=tol)

pulseSquare = pulseShaping(phiParam, pulseParam, ti, tf)

# pulseSquare.drivingStepUpPulse()
# pulseSquare.responseStepUp()

# pulseSquare.drivingStepDownPulse()
# pulseSquare.responseStepDown()

pulseSquare.drivingSquarePulse()
pulseSquare.responseSquare()
popTimeSquare = pulseSquare.getPopulationTime()
depopTimeSquare = pulseSquare.getDepopulationTime()

pulseClear = pulseShaping(phiParam, pulseParam, ti, tf, upParam, downParam)

# pulseClear.drivingClearPartPulse('up')
# pulseClear.responseClearPart('up')

# pulseClear.drivingClearPartPulse('down')
# pulseClear.responseClearPart('down')

pulseClear.drivingClearPartPulse('all')
pulseClear.responseClearPart('all')
popTimeClear = pulseClear.getPopulationTime()
depopTimeClear = pulseClear.getDepopulationTime()

fig = plt.figure(num=1, figsize=[18,8])
ax = fig.add_subplot(1, 1, 1)
ax.grid(b=True)

if pulseSquare is not None:
    ax.plot(pulseSquare.timeVec, pulseSquare.drivingPulse, ls='--', color='grey', label='Square pulse')
    ax.plot(pulseSquare.timeVec, pulseSquare.responsePulse, ls='--', color='green', label='Response to a square pulse')

if pulseClear is not None:
    ax.plot(pulseClear.timeVec, pulseClear.drivingPulse, color='black', label='CLEAR pulse')
    ax.plot(pulseClear.timeVec, pulseClear.responsePulse, color='red', label='Response to a CLEAR pulse')

if popTimeSquare is not None:
    ax.axvline(x=popTimeSquare, ls='--', color='peru')
if depopTimeSquare is not None:
    ax.axvline(x=depopTimeSquare, ls='--', color='peru')

if popTimeClear is not None:
    ax.axvline(x=popTimeClear, color='peru')
if depopTimeClear is not None:
    ax.axvline(x=depopTimeClear, color='peru')

ax.set_xlabel('Time [time unit]')
ax.set_ylabel('Amplitude of signal')
ax.set_title('Comparison between populating and depopulating time for square and CLEAR pulse')
ax.legend()
plt.show()


