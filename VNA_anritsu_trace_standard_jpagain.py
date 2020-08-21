# -*- coding: utf-8 -*-
import qt
import numpy as np
import shutil

filename = 'VNA_trace_standard' #Filename that will be used to name the data folder. Use script filename.
file_dir = 'StandardScripts/' #File directory. It has the form working_directory\file_dir\filename

###########################################################
#
#
# Sweep and measurement parameters
#
#
###########################################################

#Number of points (int) [Sample]
points = 1001

#StartFrequency [Hz]
startfrequency =6.5e9


#StopFrequency [Hz]
stopfrequency =7.5e9

#Power [dBm]
detector_power = -30.# dBm

#Detector bandwidth [Hz]
detector_bandwidth = 100. # Hz

#Averages [Sample]
averages = 1

#Sweeps
# sweeps = averages

#Sparameter
Sparam='S21'

###############################################
#
#         Coil Current Parameters
#
##############################################

#Reference coil current
# c_ref = 0e-3


################################################
#
#            Gate Voltage Parameters
#
################################################


#Reference gate voltage
# Vg_ref = 0

###########################################################
#
#
#               Devices
#
#
###########################################################


vna = qt.instruments.get('VNA')
# current_source = qt.instruments.get('hp3245a')
# Vg_source = qt.instruments.get('hp3245a_2')
# probe_src = qt.instruments.get('RFsrc')

############################################################
#
#
#       SMB100A
#
#
############################################################

# probe_src.set_status('off')

###########################################################
#
#
#               HP3245A_2 (Voltage source)
#
#
###########################################################

# Vg_source.set_mode('dcv')
# Vg_source.set_channel('A')
# Vg_source.set_resolution('high')

###########################################################
#
#
#               vna (Vector Network Analyzer)
#
#
###########################################################

vna.set_sweeptype('LIN')     #Linear frequency sweep
vna.set_trigger('IMM')       #Immediate trigger. The vna measures the second trace immediately after the first one.
# Free run in Trigger in front panel

vna.set_startfrequency(startfrequency)
vna.set_stopfrequency(stopfrequency)
vna.set_averages(averages)
vna.set_averagestatus('on')
# vna.set_sweeps(sweeps)
vna.set_port1_power(detector_power)
vna.set_port2_power(detector_power)
vna.set_measBW(detector_bandwidth)
vna.set_points(points)
vna.create_traces(('1',),(Sparam,))
vna.set_status('on')

###########################################################
#
#
#                      HP3245A
#
#
###########################################################
#
# current_source.set_mode('dci')
# current_source.set_channel('A')
# current_source.set_resolution('high')


###########################################################
#
#
#               Experiment
#
#
###########################################################

qt.mstart()

#We create the object data
data = qt.Data(name=str(filename)) #Do not change filename. Include comments as 'filename + str(comments)'

#We add the coordinate of what is swept

data.add_coordinate('Frequency [GHz]', units='GHz')

#We add the value that will be read out

data.add_value('S_21 [dB]', units='dB')
data.add_value('Phase [deg]', units='deg')
data.add_value('S_21 normed [dB]', units='dB')
data.add_value('Phase normed[deg]', units='deg')

#We live plot what we record

plot2d_1= qt.Plot2D(data, name='trace', coorddim=0, valdim=1, maxtraces=2)
plot2d_2= qt.Plot2D(data, name='trace_phase', coorddim=0, valdim=2, maxtraces=2)
plot2d_3= qt.Plot2D(data, name='trace normed', coorddim=0, valdim=3, maxtraces=2)
plot2d_4= qt.Plot2D(data, name='trace_phase normed', coorddim=0, valdim=4, maxtraces=2)



#We create a file corresponding to the object data
data.create_file()



try:

    vna.averageclear()
    # print 'here'
    vna.measure()

    qt.msleep(0.1)
    print 'here'
    ampdB, phase = vna.get_traces(('1',))[0]
    print 'here'
    qt.msleep(0.1)


    freqs=linspace(startfrequency,stopfrequency,points)/1.0e9
    data.add_data_point(freqs,ampdB,phase, ampdB-a0dB, phase - phase0)

except Exception as error:
    print str(error)
finally:
    qt.msleep()
    data.close_file()

    plot2d_1.save_png()
    plot2d_2.save_png()
    plot2d_3.save_png()
    plot2d_4.save_png()

    vna.set_status('off')
    qt.mend()
