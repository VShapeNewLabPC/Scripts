# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:48:39 2020

@author: TimGU
"""

import gc
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from lib_QND_analysis import QNDanalyse

###############################################
######### DATA FILE ###########################
###############################################

## root path
QND_folder_base = 'D:\\data_12_01_19\\Documents\\Stage_Neel\\All_data\\Seria_III\\freq_read=7.021GHz\\treadout=10\\'

## path to the folder containing the folder "Data"
QND_folder_R13_t10 = QND_folder_base + 'QND_R13_tr10\\085759___QNDness_meas_R13_tr10\\'
QND_folder_R14_t10 = QND_folder_base + 'QND_R14_tr10\\094818___QNDness_meas_R14_tr10\\'
QND_folder_R15_t10 = QND_folder_base + 'QND_R15_tr10\\101800___QNDness_meas_R15_tr10\\'
QND_folder_R16_t10 = QND_folder_base + 'QND_R16_tr10\\181444___QNDness_meas_R16_tr10\\'
QND_folder_R17_t10 = QND_folder_base + 'QND_R17_tr10\\184630___QNDness_meas_R17_tr10\\'
QND_folder_R18_t10 = QND_folder_base + 'QND_R18_tr10\\191813___QNDness_meas_R18_tr10\\'
QND_folder_R19_t10 = QND_folder_base + 'QND_R19_tr10\\195058___QNDness_meas_R19_tr10\\'
QND_folder_R20_t10 = QND_folder_base + 'QND_R20_tr10\\104528___QNDness_meas_R20_tr10\\'
QND_folder_R21_t10 = QND_folder_base + 'QND_R21_tr10\\111232___QNDness_meas_R21_tr10\\'
QND_folder_R22_t10 = QND_folder_base + 'QND_R22_tr10\\113941___QNDness_meas_R22_tr10\\'
QND_folder_R23_t10 = QND_folder_base + 'QND_R23_tr10\\120645___QNDness_meas_R23_tr10\\'
QND_folder_R24_t10 = QND_folder_base + 'QND_R24_tr10\\123348___QNDness_meas_R24_tr10\\'
QND_folder_R25_t10 = QND_folder_base + 'QND_R25_tr10\\130052___QNDness_meas_R25_tr10\\'


###############################################
######### PARAMETERS ##########################
###############################################

fit = True 
with_pio2 = False
list_QND_folder = [QND_folder_R13_t10, QND_folder_R16_t10, QND_folder_R25_t10]
list_rudat = [13, 16, 25] #dB
tread = 10 #us


# a dictionnary containing the limite of Naver to analyse
dictNaver = {13: {'all': False, 'start': 10, 'stop': 300, 'except':[]},
             14: {'all': False, 'start': 10, 'stop': None, 'except':[]},
             15: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             16: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             17: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             18: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             19: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             20: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             21: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             22: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             23: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             24: {'all': False, 'start': 10, 'stop': 1200, 'except':[]},
             25: {'all': False, 'start': 10, 'stop': 1200, 'except':[]}}

kwargs_compute = {}
## definition of the double thresold
kwargs_compute['definition'] = 'Me' #'Me' or 'Sl'
## definition of eg_center voltages
kwargs_compute['eg_center'] = 'fit' #'fit' or 'mean'
## definition of threshold to use to create the ge_sequence
kwargs_compute['threshold'] = 'double' #'double' or 'fit or 'mean'
kwargs_compute['security'] = False #True: cleaning of one point jumps
kwargs_compute['dataToPlot'] = ['rate_censor', 'rate_proba']#, 'rate_stat'] #can include 'rate', 'rate_censor', 'mean' and 'mean_censor'

dictNaverFitDecay = {13: {'start': 30 , 'stop': 200 , 'uncertainty': None},
                     14: {'start': 40 , 'stop': 260 , 'uncertainty': None},
                     15: {'start': 35 , 'stop': 700 , 'uncertainty': None},
                     16: {'start': 50 , 'stop': 800 , 'uncertainty': None},
                     17: {'start': 80 , 'stop': 420 , 'uncertainty': None},
                     18: {'start': 100 , 'stop': None , 'uncertainty': None},
                     19: {'start': 110 , 'stop': None , 'uncertainty': None},
                     20: {'start': 100 , 'stop': 400 , 'uncertainty': None},
                     21: {'start': 100 , 'stop': None , 'uncertainty': None},
                     22: {'start': 120 , 'stop': 450 , 'uncertainty': None},
                     23: {'start': 200, 'stop': 700 , 'uncertainty': None},
                     24: {'start': 280, 'stop': 1200, 'uncertainty': None},
                     25: {'start': 300, 'stop': None, 'uncertainty': None}}
dictNaverFitExcitation = {13: {'start': 30 , 'stop': 200 , 'uncertainty': None},
                          14: {'start': 15 , 'stop': None, 'uncertainty': None},
                          15: {'start': 25 , 'stop': None, 'uncertainty': None},
                          16: {'start': 50 , 'stop': 800 , 'uncertainty': None},
                          17: {'start': 80 , 'stop': 420 , 'uncertainty': None},
                          18: {'start': 150, 'stop': 700 , 'uncertainty': None},
                          19: {'start': 220, 'stop': 900 , 'uncertainty': None},
                          20: {'start': 230, 'stop': 430 , 'uncertainty': None},
                          21: {'start': 300, 'stop': None , 'uncertainty': None},
                          22: {'start': 440, 'stop': None, 'uncertainty': None},
                          23: {'start': 250, 'stop': None, 'uncertainty': None},
                          24: {'start': 270, 'stop': 900 , 'uncertainty': None},
                          25: {'start': 330, 'stop': 800, 'uncertainty': None}}
dictNaverFit = {'decay': dictNaverFitDecay, 'excitation': dictNaverFitExcitation}      
         
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
print('Size of objects at the start : {}'.format(sizeof_start))

plt.close('all')
start = time.time()
     
anaQND = QNDanalyse(list_QND_folder, list_rudat, dictNaver, with_pio2, **kwargs_compute)
anaQND.loadQNDresult()

anaQND.setListNaver_fromMaxQe()
savepath = anaQND.getSavePath(list_QND_folder, QND_folder_base)

if fit:
    anaQND.fitAllJumpRate(dictNaverFit) 
    anaQND.plotJumpRateFitVsRudat(tread, dictNaverFit, RudatToNbPhoton=(22,1.6))#, savepath=savepath)



print('')
print("###################")
print("Start plot vs naver")

if fit:
    start = time.time()
    for rudat in anaQND.listRudat:
        sys.stdout.write('\rPlot JumpRate vs naver fit : rudat = {}dB'.format(rudat))
        anaQND.plotJumpRateVsNaver(rudat, fit=True, savepath=savepath+'JumpRate_vs_naver_fit\\', savename='Jump_rate_vs_naver_R{}_tr10.png'.format(rudat))
    
    end = time.time()
    sys.stdout.write('\rPlot JumpRate vs naver fit ... done - {:.3f}s\n'.format(end-start))

##############################

start = time.time()
for rudat in anaQND.listRudat:
    sys.stdout.write('\rPlot JumpRate vs naver raw : rudat = {}dB'.format(rudat))
    anaQND.plotJumpRateVsNaver(rudat, fit=False, savepath=savepath+'JumpRate_vs_naver_raw\\', savename='Jump_rate_vs_naver_R{}_tr10.png'.format(rudat))

end = time.time()
sys.stdout.write('\rPlot JumpRate vs naver raw ... done - {:.3f}s\n'.format(end-start))
   
##############################

start = time.time()
for rudat in anaQND.listRudat:
    sys.stdout.write('\rPlot QNDness vs naver : rudat = {}dB'.format(rudat))
    anaQND.plotQNDnessVsNaver(rudat, savepath=savepath+'QNDness_vs_naver\\', savename='QNDness_vs_naver_R{}_tr10.png'.format(rudat))

end = time.time()
sys.stdout.write('\rPlot QNDness vs naver ... done - {:.3f}s\n'.format(end-start))



############################


print('')
print("###################")
print("Start plot vs rudat")
start = time.time()
for Naver in anaQND.listAllNaver:
    sys.stdout.write('\rPlot JumpRate vs rudat : n_aver = {}ns'.format(Naver))
    anaQND.plotJumpRateVsRudat(Naver, savepath=savepath+'JumpRate_vs_rudat\\', savename='Jump_rate_vs_rudat_naver{}_tr10.png'.format(Naver))

end = time.time()
sys.stdout.write('\rPlot JumpRate vs rudat ... done - {:.3f}s\n'.format(end-start))

##############################

start = time.time()
for Naver in anaQND.listAllNaver:
    sys.stdout.write('\rPlot QNDness vs rudat : n_aver = {}ns'.format(Naver))
    anaQND.plotQNDnessVsRudat(Naver, savepath=savepath+'QNDness_vs_rudat\\', savename='QNDness_vs_rudat_naver{}_tr10.png'.format(Naver))

end = time.time()
sys.stdout.write('\rPlot QNDness vs rudat ... done - {:.3f}s\n'.format(end-start))

    
print('Size of objects at the start : {}'.format(sizeof_start))
print('size of objects at the end: {}'.format(sys.getsizeof(gc.get_objects())))
   
end = time.time()
print('\nTotal duration = {}s'.format(get_duration(end-start)))
print('end')    