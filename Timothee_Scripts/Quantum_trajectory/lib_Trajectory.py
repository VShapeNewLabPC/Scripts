# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:49:53 2020

@author: TimGU
"""

import gc

import numpy as np
import matplotlib.pyplot as plt


STATE_G = 0
STATE_E = 1

MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def vectorCompression(vector, indStart, indEnd, n):
    '''
    Decrease the length of given vector by averaging neighbor values
    Takes only the part between [ind_start, ind_end]
    n2avr - number of point to average
    Example (witn n2avr=2):
    [0,2,3,5] -> [1, 4]
    
    Example (witn n2avr=3):
    [1,2,3,4,5,6,7,4] -> [2.  5.  5.5]
    '''
    if n == 0:
        return vector
    
    if n == 1:
        return vector[indStart:indEnd]
    
    vector = vector[indStart:indEnd]
    
    if len(vector)%n != 0:
        padVector = np.pad(vector.astype(float), (0, n - vector.size%n), 
                           mode='constant', constant_values=np.NaN).reshape(-1, n)
        meanVector = np.nanmean(padVector, axis=1)
    else:
        padVector = vector.reshape(-1, n)
        meanVector = np.mean(padVector, axis=1)
    
    return meanVector


 

############################################################################
###################   CLASS TRAJECT   ######################################
############################################################################

class Traject(object):
    '''
    Class for containing and characterizing a single trajectorie

    Parameters
    ----------
    type_trj : string
        Either 'no', 'pi' or 'pio2'.
    args_data : list
        Contains [t_vec, I_vec, Q_vec].
    args_loading : list
        Contains [n_aver].
    '''
    
    type_trj = '' # 'pi', 'no' or 'pio2'

    ## raw data ##
    time = None
    re = None
    im = None

    ## normalized, compressed data ##
    time_comp = None
    x_comp = None

    ## digital list (g or e state) ##
    ge_seq = None
    n_jump = 0
    n_not_jump = 0

    p_gg = 0
    p_ee = 0
    p_ge = 0
    p_eg = 0
    
    ## thesholds ##
    threshold = None
    threshold_low = None
    threshold_high = None
    
    g_center = None
    e_center = None
    g_size = None
    e_size = None

    n_aver = None
    

    ####______Methods______________####
    def __init__(self, type_trj, args_data, args_loading):
        '''
        Initialisation method of the class Traject

        Parameters
        ----------
        type_trj : string
            Either 'no', 'pi' or 'pio2'.
        args_data : list
            Contains [t_vec, I_vec, Q_vec].
        args_loading : list
            Contains [n_aver].

        Returns
        -------
        None.

        '''

        t_vec, I_vec, Q_vec = args_data[0], args_data[1], args_data[2]
        
        self.time = t_vec
        if I_vec is not None:
            self.re = I_vec*1000 #conversion from [V] to [mV]
        if Q_vec is not None:
            self.im = Q_vec*1000 #conversion from [V] to [mV]

        self.type_trj = type_trj

        self.n_aver = args_loading[0]


    ################################
    ###   SET TRAJ  ################
    def setNormalizedData(self, **kwargs):
        '''
        Rotate and compressed data

        Parameters
        ----------
        **kwargs : dict 
            Contains keys 'n_aver', 'ind_start' and 'ind_stop'

        Returns
        -------
        None.

        '''
        self.n_aver = kwargs['n_aver'] 
        ind_start = kwargs['ind_start'] 
        ind_stop = kwargs['ind_stop']
        
        
        if self.time is not None:
            self.time_comp =  vectorCompression(self.time, ind_start, ind_stop, self.n_aver)
        if self.re is not None and self.im is not None:
            covMat = np.cov(self.re, self.im)
            eigenValue, eigenVector = np.linalg.eig(covMat)
            X = -eigenVector[0,0]*self.re + eigenVector[0,1]*self.im
            # Y =  eigenVector[1,0]*self.re - eigenVector[1,1]*self.im
            
            self.x_comp = vectorCompression(X, ind_start, ind_stop, self.n_aver)
     
        
    def getState_SingleThreshold(self, val, threshold, sign_ge):
        '''
        Gives the state of the points using only one threshold

        Parameters
        ----------
        val : float
            the voltage value of the point.
        threshold : float
            the threshold value to use.
        sign_ge : int
            describe whether |e> state is above |g> state (positive value) 
            or not (negative value)

        Returns
        -------
        STATE
            Either STATE_E or STATE_G.

        '''
        if sign_ge > 0:
            if val > threshold:
                return STATE_E
            else:
                return STATE_G
        else:
            if val > threshold:
                return STATE_G
            else:
                return STATE_E
        
    
    def setStateSeq_singleThreshold(self, **kwargs):
        '''
        Create a sequence of g or e states using one threshold

        Parameters
        ----------
        **kwargs : dict 
            Contains keys 'sign_ge', 'threshold', 'security'

        Returns
        -------
        None.

        '''

        sign_ge = kwargs['sign_ge']
        threshold = kwargs['threshold']
        security = kwargs['security']

        if self.x_comp is None:
            print('cant set ge_sequence -- x_comp does not exist')
            return False
        if threshold is not None:
            self.threshold = threshold

        ###-------------###
        seq_state = np.zeros_like(self.x_comp)
        for i in range(len(self.x_comp)):
            seq_state[i] = self.getState_SingleThreshold(self.x_comp[i],  threshold, sign_ge)
            
        if security:
            #Clean_single_points
            for i in range(len(seq_state)-1):
                if i == 0:
                    continue
                    # we dont check first and last point here
                if seq_state[i+1] == seq_state[i-1]:
                    if seq_state[i] != seq_state[i+1]:
                        seq_state[i] = seq_state[i+1]
        ## now work on edge points:
        ## sure?
        seq_state[0] = seq_state[1]
        seq_state[-1] = seq_state[-2]

        ###-------------###
        ## IMPROVE: make killing more than one point of fluctuation

        self.ge_seq = seq_state
    
    
    def getState_doubleThreshold(self, val, previous_state, threshold_low, threshold_high, sign_ge):
        '''
        Gives the state of the points using double threshold technic

        Parameters
        ----------
        val : float
            The voltage value of the point.
        previous_state : int
            Either STATE_E or STATE_G.
        threshold_low : float
            The value of the low threshold.
        threshold_high : float
            The value of the high threshold.
        sign_ge : int
            describe whether |e> state is above |g> state (positive value) 
            or not (negative value)

        Returns
        -------
        STATE
            Either STATE_E or STATE_G.

        '''
        if sign_ge > 0:
            if val > threshold_high:
                return STATE_E
            elif val < threshold_low:
                return STATE_G
            else:
                return previous_state
        else:
            if val > threshold_high:
                return STATE_G
            elif val < threshold_low:
                return STATE_E
            else:
                return previous_state


    def setStateSeq_doubleThreshold(self, **kwargs):
        '''
        Create a sequence of g or e states using double threshold technic

        Parameters
        ----------
        **kwargs : dict 
            Contains keys 'sign_ge', 'threshold_fit', 'threshold_low', 
            'threshold_high' and 'security'

        Returns
        -------
        None.

        '''
        sign_ge = kwargs['sign_ge']
        threshold_fit = kwargs['threshold_fit']
        threshold_low = kwargs['threshold_low']
        threshold_high = kwargs['threshold_high']
        security = kwargs['security']
        
        seq_state = np.zeros_like(self.x_comp)
            
        if self.type_trj == 'pi':
            seq_state[0] = STATE_E
        elif self.type_trj == 'no':
            seq_state[0] = STATE_G
        else:
            seq_state[0] = self.getState_SingleThreshold(self.x_comp[0],  threshold_fit, sign_ge)
            
        for i in range(1, len(self.x_comp)):
            seq_state[i] = self.getState_doubleThreshold(self.x_comp[i], seq_state[i-1], 
                                                threshold_low, threshold_high, sign_ge)
            
        if security:
            #Clean_single_points
            for i in range(len(seq_state)-1):
                if i == 0:
                    continue
                    # we dont check first and last point here
                if seq_state[i+1] == seq_state[i-1]:
                    if seq_state[i] != seq_state[i+1]:
                        seq_state[i] = seq_state[i+1]
            
        self.ge_seq = seq_state
          

    def setJumpCounter(self, **kwargs):
        '''
        Count the number of jumps, compute Pgg, Pee, Pge and Peg

        Parameters
        ----------
        **kwargs : dict
            Not used.

        Returns
        -------
        None.

        '''
        if self.ge_seq is None:
            print('there is no ge_seq to count jumps')
            return False

        self.n_jump = 0
        self.n_not_jump = 0

        ## just try to check remys way !V bullshit..
        self.p_gg = 0
        self.p_ee = 0
        self.p_ge = 0
        self.p_eg = 0
        for i in range( len(self.ge_seq)-1 ):
            if   (self.ge_seq[i] == STATE_G) and (self.ge_seq[i+1] == STATE_G):
                self.p_gg +=1
                self.n_not_jump +=1
            elif (self.ge_seq[i] == STATE_E) and (self.ge_seq[i+1] == STATE_E):
                self.p_ee +=1
                self.n_not_jump +=1
            elif (self.ge_seq[i] == STATE_G) and (self.ge_seq[i+1] == STATE_E):
                self.p_ge +=1
                self.n_jump +=1
            elif (self.ge_seq[i] == STATE_E) and (self.ge_seq[i+1] == STATE_G):
                self.p_eg +=1
                self.n_jump +=1
            else:
                print('WTF?? jump_counter()')
    
    
    def plotter(self, savepath='', savename=''):
        '''
        Plots compressed and rotated(in right axes) trajectorie 
        Optionnal plot ge-sequence, e and g centers and size, thresholds

        Parameters
        ----------
        savepath : str, optional
            The path to the saving folder. The default is ''.
        savename : TYPE, optional
            The name . The default is ''.

        Returns
        -------
        bool
            If the plot has been made return True, else return False

        '''

        if self.x_comp is None:
            print('Plot of trajectory impossible: no data to plot')
            return False
        
        fig = plt.figure(figsize=[10,7])
        ax1 = fig.add_subplot(1, 1, 1)

        ax1.plot(self.time_comp, self.x_comp, color='r')
        
        if self.threshold is not None:
            ax1.axhline(y=self.threshold, color='peru')
        
        if self.threshold_low is not None:
            ax1.axhline(y=self.threshold_low, color='black')
            
        if self.threshold_high is not None:
            ax1.axhline(y=self.threshold_high, color='black')

        if self.g_center is not None:
            ax1.axhline(y=self.g_center, color='b')

            if self.g_size is not None:
                ax1.fill_between(self.time_comp,  self.g_center-self.g_size/2,
                                                  self.g_center+self.g_size/2, color='b', alpha=0.4)
        if self.e_center is not None:
            ax1.axhline(y=self.e_center, color='r')

            if self.e_size is not None:
                ax1.fill_between(self.time_comp,  self.e_center-self.e_size/2,
                                                  self.e_center+self.e_size/2, color='r', alpha=0.4)

        if self.ge_seq is not None:
            ax2 = ax1.twinx()
            ax2.plot(self.time_comp, self.ge_seq, color='grey', label='g-e-state', 
                     alpha=1.0, lw=3)
            ax2.set_yticks([0, 1])
            ax2.set_yticklabels(['g state', 'e state'])

        ax1.set_xlabel('Time [ns]')
        ax1.set_ylabel('Voltage [mV]')

        fig.tight_layout()

        if savepath != '':
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)

            plt.savefig(savepath+savename)
            plt.cla()
            plt.clf()
            plt.close(fig)
            gc.collect()
        else:
            plt.show()
        
        return True
