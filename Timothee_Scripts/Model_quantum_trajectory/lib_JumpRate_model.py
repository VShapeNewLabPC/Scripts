# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 08:54:37 2020

@author: TimGU
"""

import gc
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import gamma
from scipy.optimize import curve_fit
from scipy import integrate

from lib_QND_processingModel import QNDness

T1 = 2000 #ns
STATE_E = 1
STATE_G = 0

def get_duration(duration):
    duration = int(duration)
    s = duration%60
    m = duration//60
    if m >= 60:
        m = m%60
    h = duration//3600
    
    return '{:02d}:{:02d}:{:02d}'.format(h, m, s)


def take_dt_str(date=True, hour=True, slash=True):
    import datetime
    s1 = ''
    s = str(datetime.datetime.now())

    if date:
        s1 =s1 +s[2]+s[3]+s[5]+s[6]+s[8]+s[9]
    if slash:
        s1 = s1 + '_'
    if time:
        s1 = s1 + s[11]+s[12]+s[14]+s[15]

    return s1



def gammaDistribution(t, amp, rate, alpha):
    return amp * (rate**alpha) * (t**(alpha-1)) * np.exp(-rate*t) / gamma(alpha)

def gammaDistributionLabel(amp, rate, alpha):
    label = r'$%.3f$ '%(amp)
    label += r'$%.3f^{%.3f}$'%(rate, alpha)
    label += r'$t^{%.3f - 1}$'%(alpha)
    label += r'$exp({%.3f t})$'%(rate)
    label += r'$/\Gamma({%.3f})$'%(alpha)
    return label


######################################################################
        
def computeHistogram(data, nbins=None, normalized=True, verify=None):
    if nbins is None:
        bins = 'auto'
        if verify is not None:
            n = len(data)
            iqr = np.subtract(*np.percentile(data, [75, 25]))
            h = 2*iqr/n**(1/3)
            print(h, verify)
            if h<verify:
                rangeData = np.max(data) - np.min(data)
                bins = int(np.round(np.ceil(rangeData / verify)))
    else:
        bins = nbins
   
    data_hist, bin_edges = np.histogram(data, bins=bins)
    binscenters_hist = np.array([0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(len(bin_edges)-1)])
    if normalized:
        data_hist = data_hist / len(data)
    
    return data_hist, binscenters_hist, bin_edges


######################################################################

class modelJumpRate(QNDness):
        
    listNaver = None
    dictNbTraject = None
    
    durationTraject = None
    deltaStart = None
    deltaStop = None
    
    decay_rate = None
    excitation_rate = None
    
    coeffGaussian_e = None
    coeffGaussian_g = None
    coeffGaussian_delta = None
    
    def __init__(self, **kwargs):
        print('Initialisation of model')
        
        self.dictNbTraject = {}
        self.dictNbTraject['no'] = kwargs['nb_trj_no']
        self.dictNbTraject['pi'] = kwargs['nb_trj_pi']
        
        self.durationTraject = kwargs['durationTraject'] #ns
        self.deltaStart = kwargs['deltaStart'] #ns
        self.deltaStop = kwargs['deltaStop'] #ns
        
        
        self.decay_rate = kwargs['decay_rate']/1000 #conversion from /us to /ns
        self.excitation_rate = kwargs['excitation_rate']/1000 #conversion from /us to /ns
        
        self.coeffGaussian_e = kwargs['coeffGaussian_e']
        self.coeffGaussian_g = kwargs['coeffGaussian_g']
        self.coeffGaussian_delta = kwargs['coeffGaussian_delta']
        
        self.listNaver = kwargs['listNaver']
        
        savepath = kwargs['savepath']
        saveFolder = 'ER{}_DR{}_'.format(self.excitation_rate*1000, self.decay_rate*1000) + 'Traject_nb{}_Dt{}us_N{}\\'.format(self.dictNbTraject['no'], 
            self.durationTraject/1000, self.coeffGaussian_e[1])
        QNDness.__init__(self, savepath+saveFolder, False, self.durationTraject/1000, 0, None)

    
    def createRe_Traject(self, type_trj):
        Re = []
        time_trj = 0
        
        if type_trj == 'pi':
            state_init = STATE_E
            next_state = STATE_G
        else:
            state_init = STATE_G
            next_state = STATE_E
        
        for _ in range(self.deltaStart):
            time_trj += 1
            voltage = np.random.normal(*self.coeffGaussian_delta)
            Re.append(voltage/1000.0) #conversion from mV to V
        
        while time_trj < self.durationTraject:
            if state_init == STATE_G:
                stateDuration = np.random.exponential(scale=1/self.excitation_rate)
            else:
                stateDuration = np.random.exponential(scale=1/self.decay_rate)
                
            if time_trj+stateDuration >= self.durationTraject:
                stateDuration = self.durationTraject - time_trj
            
            dwell_time = 0
            while dwell_time < stateDuration:
                if state_init == STATE_G:
                    voltage = np.random.normal(*self.coeffGaussian_g)
                else:
                    voltage = np.random.normal(*self.coeffGaussian_e)
                    
                Re.append(voltage/1000.0) #conversion from mV to V
                dwell_time += 1
                time_trj += 1
                
            state_init, next_state = next_state, state_init
        
        for _ in range(self.deltaStop):
            time_trj += 1
            voltage = np.random.normal(*self.coeffGaussian_delta)
            Re.append(voltage/1000.0) #conversion from mV to V
        
        return np.array(Re)
        

    def createAllTraject(self):
        self.loading_fromRawData = True
        self.loading_fromProcessedData = False
        
        time_vec = np.arange(0, self.durationTraject+self.deltaStop, 1) #ns
        Im = np.zeros(np.shape(time_vec))
        for type_trj in ['pi', 'no']:
            
            start = time.time()
            n_max = self.dictNbTraject[type_trj]
            for n in range(n_max):
                sys.stdout.write('\r{} : creation of raw trajctories {}/{}'.format(type_trj, n+1, n_max))
                
                Re = self.createRe_Traject(type_trj)
                args_data = [time_vec, Re, Im]
                args_loading = [0]
                self.addTraject_fromData(type_trj, args_data, args_loading)
            
            end = time.time()
            sys.stdout.write('\r{} : creation of raw trajctories     ... done - {:.3f}s\n'.format(type_trj, end-start))
      
        
    def saveAllTraject(self):
        savepath = self.dataFolder[0] + 'Data\\'
        import os
        if not os.path.exists(savepath):
            os.makedirs(savepath)
     
        for type_trj in ['pi', 'no']:
            start = time.time()
            sys.stdout.write('{} : saving of raw trajctories       ...'.format(type_trj))
            
            header = 'time [ns] \t I_{}(t) [V] \t Q_{}(t) [V]'.format(type_trj, type_trj)
            
            ntrj, nfile = 1, 0
            listAllRe = np.array([])
            listAllIm = np.array([])
            listAllTime = np.array([])
            for trj in self.dictAllTraject[type_trj]:

                listAllRe = np.concatenate((listAllRe, trj.re/1000))
                listAllIm = np.concatenate((listAllIm, trj.im/1000))
                listAllTime = np.concatenate((listAllTime, trj.time))
                    
                if ntrj == 100:
                    nfile += 1
                    ntrj = 1
                    
                    savename = 'traj_{}_{}_{}.dat'.format(type_trj, nfile, take_dt_str(date=False, hour=True, slash=False))
                    mat_data = np.column_stack((listAllTime, listAllRe, listAllIm))
                    np.savetxt(savepath+savename, mat_data, header=header)
                    
                    listAllRe = np.array([])
                    listAllIm = np.array([])
                    listAllTime = np.array([])
                else:
                    ntrj += 1
             
            if ntrj < 100:
                nfile += 1
                
                savename = 'traj_{}_{}_{}.dat'.format(type_trj, nfile, take_dt_str(date=False, hour=True, slash=False))
                mat_data = np.column_stack((listAllTime, listAllRe, listAllIm))
                np.savetxt(savepath+savename, mat_data, header=header)
                
            end = time.time()
            sys.stdout.write('\r{} : saving of raw trajctories       ... done - {:.3f}s\n'.format(type_trj, end-start))  
    
        
    
    def proceedSimulatedData(self, list_pplot, **kwargs_compute):
        save_result = kwargs_compute['save_result']
        
        start = time.time()
        for n_aver, pplot in zip(self.listNaver, list_pplot):
            print('\n############')
            print('n_aver = {}, excitation rate = {}/us, decay rate= {}/us'.format(n_aver, self.excitation_rate*1000, self.decay_rate*1000))
            start_loop = time.time()
    
            self.loadingProcedure(n_aver)
              
            self.compute_all(**kwargs_compute)
            self.save_all(percent_plot=pplot, **kwargs_compute)
            
            end_loop = time.time()
            print(' ')
            print('duration of one loop of the process = {}s'.format(get_duration(end_loop-start_loop)))
            print('Garbage cleaner in loop: {}+{}+{}'.format(gc.collect(), gc.collect(), gc.collect()))
                
        end = time.time()
        
        print(' ')
        print('duration of all the process = {}s'.format(get_duration(end-start)))
        print('Garbage cleaner end of process: {}+{}+{}'.format(gc.collect(), gc.collect(), gc.collect()))
            
                
        
   
    
      
        


        
        
    