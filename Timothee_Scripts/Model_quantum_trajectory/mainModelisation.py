# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:48:39 2020

@author: TimGU
"""

import gc
import sys
import time

import numpy as np

from lib_JumpRate_model import modelJumpRate


def get_duration(duration):
    duration = int(duration)
    s = duration%60
    m = duration//60
    if m >= 60:
        m = m%60
    h = duration//3600
    
    return '{:02d}:{:02d}:{:02d}'.format(h, m, s)


sizeof_start = sys.getsizeof(gc.get_objects())
start = time.time()


###############################################
######### PARAMETERS ##########################
###############################################


kwargs_compute = {}
# definition of the double thresold
kwargs_compute['definition'] = 'Me' #'Me' or 'Sl'
# definition of eg_center voltages
kwargs_compute['eg_center'] = 'fit' #'fit' or 'mean'
#definition of threshold to use to create the ge_sequence
kwargs_compute['threshold'] = 'double' #'double' or 'fit or 'mean'
kwargs_compute['security'] = False #True: cleaning of one point jumps
kwargs_compute['size'] = True #True: give size of ge_state to trajectory for plotter
kwargs_compute['save_result'] = True #True: save all figures
kwargs_compute['saveResult'] = True #True: save all figures
#True: save only processed data, for first processng of data, can be good to only
#save processed data, avoid time consumming problems
kwargs_compute['saveOnlyProcessed'] = False 


list_naver = np.concatenate((np.arange(0, 50, 5), np.arange(50, 451, 10), np.arange(500, 2001, 100)))  

## list of the averaging number for which plot of individual quantum trajectory
## will be perform
list_NaverPlot = [0, 50, 100] # I want that for n_aver=0, 50 and 100 plot of individual trajectory is perform

## the number of plot to do in percent of the total number of trajectory
nbPlot_perNaver = 10 #%

# set the list of plot percent
list_pplot = np.zeros(np.shape(listNaver))
for Naver in list_NaverPlot:
    list_pplot[list_naver == Naver] = nbPlot_perNaver


kwargs_loading = {}
kwargs_loading['listNaver'] = list_naver
kwargs_loading['savepath'] = 'D:\\data_12_01_19\\Documents\\Stage_Neel\\Modelisation\\'

#before the pi-pulse if applied, measurement is started
kwargs_loading['deltaStart'] = 150 #ns
#after the measurement is stopped, continue to acquire data
kwargs_loading['deltaStop'] = 230 #ns

kwargs_loading['coeffGaussian_e'] = (10.0, 10.0) #mean and std
kwargs_loading['coeffGaussian_g'] = (-10.0, 11.0) #mean and std
kwargs_loading['coeffGaussian_delta'] = (0.0, 10.0) #mean and std

kwargs_loading['decay_rate'] = 2 #/us
kwargs_loading['excitation_rate'] = 0.15 #/us

kwargs_loading['nb_trj_no'], kwargs_loading['nb_trj_pi'] = 1000, 1000
kwargs_loading['durationTraject'] = 10000 #ns


#######################################################################
######################---------------------------######################
######################___________MAIN____________######################
######################---------------------------######################
#######################################################################

model = modelJumpRate(**kwargs_loading)
model.createAllTraject()
model.proceedSimulatedData(list_pplot, **kwargs_compute)


print('Size of objects at the start : {}'.format(sizeof_start))
print('size of objects at the end: {}'.format(sys.getsizeof(gc.get_objects())))
   
end = time.time()
print('\nTotal duration = {}s'.format(get_duration(end-start)))
print('end')    



    