# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:48:39 2020

@author: TimGU
"""

import gc
import sys
import time

import numpy as np

from datetime import datetime

from lib_QND_processing import QNDness


###############################################
######### DATA FILE ###########################
###############################################

## root path
# QND_folder_t10 = 'D:\\data_12_01_19\\Documents\\Stage_Neel\\All_data\\Seria_III\\freq_read=7.021GHz\\treadout=10\\'

QND_folder_t10 = 'C:\\Users\\VShapeExperiment\\Documents\\Timothee_QND_Data\\TrajectoryBeforeConfinement\\'

## path to the folder containing the folder "Data"
QND_folder_R13_t10 = QND_folder_t10 + 'QND_R13_tr10\\085759___QNDness_meas_R13_tr10\\'
QND_folder_R14_t10 = QND_folder_t10 + 'QND_R14_tr10\\094818___QNDness_meas_R14_tr10\\'
QND_folder_R15_t10 = QND_folder_t10 + 'QND_R15_tr10\\101800___QNDness_meas_R15_tr10\\'
QND_folder_R16_t10 = QND_folder_t10 + 'QND_R16_tr10\\181444___QNDness_meas_R16_tr10\\'
QND_folder_R17_t10 = QND_folder_t10 + 'QND_R17_tr10\\184630___QNDness_meas_R17_tr10\\'
QND_folder_R18_t10 = QND_folder_t10 + 'QND_R18_tr10\\191813___QNDness_meas_R18_tr10\\'
QND_folder_R19_t10 = QND_folder_t10 + 'QND_R19_tr10\\195058___QNDness_meas_R19_tr10\\'
QND_folder_R20_t10 = QND_folder_t10 + 'QND_R20_tr10\\104528___QNDness_meas_R20_tr10\\'
QND_folder_R21_t10 = QND_folder_t10 + 'QND_R21_tr10\\111232___QNDness_meas_R21_tr10\\'
QND_folder_R22_t10 = QND_folder_t10 + 'QND_R22_tr10\\113941___QNDness_meas_R22_tr10\\'
QND_folder_R23_t10 = QND_folder_t10 + 'QND_R23_tr10\\120645___QNDness_meas_R23_tr10\\'
QND_folder_R24_t10 = QND_folder_t10 + 'QND_R24_tr10\\123348___QNDness_meas_R24_tr10\\'
QND_folder_R25_t10 = QND_folder_t10 + 'QND_R25_tr10\\130052___QNDness_meas_R25_tr10\\'


###############################################
######### PARAMETERS ##########################
###############################################


with_pio2 = False
list_QND_folder = [QND_folder_R13_t10, QND_folder_R16_t10, QND_folder_R25_t10]
list_rudat = [13, 16, 25] #dB
list_tread = [10, 10, 10] #us

## if processed data already exist, list_naver can equal 'all'
## if not, or want a specific Naver, can set it manually
# list_naver = 'all'
list_naver = np.concatenate((np.arange(0, 50, 5), np.arange(50, 451, 10), np.arange(500, 2001, 100)))

## list of the averaging number for which plot of individual quantum trajectory
## will be perform
list_NaverPlot = [0, 50, 100] # I want that for n_aver=0, 50 and 100 plot of individual trajectory is perform

## the number of plot to do in percent of the total number of trajectory
nbPlot_perNaver = 10 #%


kwargs = {}
## definition of the double thresold
kwargs['definition'] = 'Me' #'Me' or 'Sl'
## definition of eg_center voltages
kwargs['eg_center'] = 'fit' #'fit' or 'mean'
## definition of threshold to use to create the ge_sequence
kwargs['threshold'] = 'double' #'double' or 'fit or 'mean'
kwargs['security'] = False #True: cleaning of one point jumps
kwargs['size'] = True #True: give size of ge_state to trajectory for plotter
kwargs['saveResult'] = False #True: save all figures
## True: save only processed data, for first processng of data, can be good to only
## save processed data, avoid time consumming problems
kwargs['saveOnlyProcessed'] = False

def get_duration(duration):
    duration = int(duration)
    s = duration%60
    m = duration//60
    if m >= 60:
        m = m%60
    h = duration//3600

    return '{:02d}:{:02d}:{:02d}'.format(h, m, s)


#######################################################################
######################---------------------------######################
######################___________MAIN____________######################
######################---------------------------######################
#######################################################################

print('Garbage cleaner start : {}'.format(gc.collect()))

sizeof_start = sys.getsizeof(gc.get_objects())
print('Size of objects in RAM at the start : {}'.format(sizeof_start))

###############################################

start_all = time.time()

# loop over the list of QND folder
for QND_folder, tread, rudat in zip(list_QND_folder, list_tread, list_rudat):
    start = time.time()

    meaQND = QNDness(QND_folder, with_pio2, tread, rudat)

    if type(list_naver) is str and list_naver == 'all':
        # set automatically the list of averaging number by looking at the folder
        # "Processed_data"
        listNaver = meaQND.getNaverList()
    else:
        listNaver = list_naver


    # set the list of plot percent
    list_pplot = np.zeros(np.shape(listNaver))
    for Naver in list_NaverPlot:
        list_pplot[listNaver == Naver] = nbPlot_perNaver

    # loop over the number of averaging
    for n_aver, pplot in zip(listNaver, list_pplot):
        print('\n############')
        print('n_aver = {}, rudat = {}dB'.format(n_aver, rudat))
        start_loop = time.time()

        #############################
        ######### MAIN PART #########
        #############################

        #===========================================================
        # LOADING PROCEDURE
        # if data loaded from raw data:
            # normalized it using n_aver
        # else:
            # try to load processed data
            # if fail:
                # load raw data and normalized it using n_aver
        meaQND.loadingProcedure(n_aver)
        #===========================================================
        #===========================================================
        # COMPUTATION PROCEDURE (in parenthesis: the class doing the computation)
        # compute statistical trajectories (MeanTraject)
        # compute single threshold using statistical trajectories (MeanTraject)
        # compute and fit by double gaussian histogram of all voltage points (BigHistogram)
        # compute double threshold using the double fit (BigHistogram)
        # compute state sequence of all trajectories (Traject)
        # compute duration list of staying in e or g states (StateDuration)
        # compute jump rates by fitting the histogram of duration (StateDuration)
        # compute voltage histogram of e and g states (StateVoltage)
        # compute probability to be in e or g states at each time (MeanProbaState)
        # compute jump rates by fitting the probability (MeanProbaState)
        # compute QNDness (QNDness)
        meaQND.compute_all(**kwargs)
        #===========================================================
        #===========================================================
        meaQND.save_all(percent_plot=pplot, **kwargs)

        #############################

        end_loop = time.time()
        print(' ')
        print('duration of one loop of the process = {}s'.format(get_duration(end_loop-start_loop)))
        print('Garbage cleaner in loop: {}+{}+{}'.format(gc.collect(), gc.collect(), gc.collect()))

    end = time.time()

    print(' ')
    print('duration of all the process = {}s'.format(get_duration(end-start)))
    print('Garbage cleaner end of process: {}+{}+{}'.format(gc.collect(), gc.collect(), gc.collect()))

    sizeof_end = sys.getsizeof(gc.get_objects())
    print('Size of objects in RAM at the start : {}'.format(sizeof_start))
    print('Size of objects in RAM at the end: {}'.format(sizeof_end))
    if sizeof_end != sizeof_start:
        print('Size of object in the RAM is not the same at the start and at the end of the processing, problem of memory leakage !!!')

end_all = time.time()
print('\nTotal duration = {}s'.format(get_duration(end_all-start_all)))
print('end at {}'.format(str(datetime.now())))

########################################################################################
