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
points = 2001

#StartFrequency [Hz]
startfrequency =6.5e9


#StopFrequency [Hz]
stopfrequency =8.5e9

#Power [dBm]
detector_power = -0.# dBm

#Detector bandwidth [Hz]
detector_bandwidth = 1000 # Hz

#Averages [Sample]
averages = 1

#Sweeps
sweeps = averages

#Sparameter
Sparam='S21'

###############################################
#
#         Coil Current Parameters
#
##############################################

#Reference coil current
c_ref =-0.0


################################################
#
#            Gate Voltage Parameters
#
################################################


#Reference gate voltage
Vg_ref = 0

###########################################################
#
#
#               Devices
#
#
###########################################################


vna = qt.instruments.get('vna')
# current_source = qt.instruments.get('hp3245a')
# Vg_source = qt.instruments.get('hp3245a_2')
#probe_src = qt.instruments.get('probe_src')

############################################################
#
#
#       SMB100A
#
#
############################################################

#probe_src.set_status('off')

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
#               ZNB (Vector Network Analyzer)
#
#
###########################################################

vna.initialize_one_tone_spectroscopy(('trace1',),(Sparam,))

vna.set_startfrequency(startfrequency)
vna.set_stopfrequency(stopfrequency)
vna.set_averages(averages)
vna.set_averagestatus('on')
vna.set_sweeps(sweeps)
vna.set_power(detector_power)
vna.set_measBW(detector_bandwidth)
vna.set_points(points)


###########################################################
#
#
#                      HP3245A
#
#
###########################################################

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

data.add_value(str(Sparam)+'[dB]', units='dB')
data.add_value('Phase [deg]', units='deg')


#We live plot what we record

plot2d_1= qt.Plot2D(data, name='trace', coorddim=0, valdim=1, maxtraces=2)
plot2d_2= qt.Plot2D(data, name='trace_phase', coorddim=0, valdim=2, maxtraces=2)




#We create a file corresponding to the object data
data.create_file()



try:
    # Vg_source.on()
    # Vg_source.set_voltage(Vg_ref)
    # print 'here'
    # current_source.on()
    # qt.msleep(2)
    # current_source.set_current(c_ref)
    # qt.msleep(2)

    vna.averageclear()
    # print 'here'
    vna.measure()
    qt.msleep(0.1)
    amp,phase=vna.get_traces( ('trace1',))[0]
    qt.msleep(0.1)
    # print 'here'

    freqs=linspace(startfrequency,stopfrequency,points)/1.0e9
    data.add_data_point(freqs,amp,phase)

except Exception as error:
    print str(error)
finally:
    qt.msleep()
    data.close_file()
    data_dir=qt.Data.get_dir(data)           #Where to save the copy of the script.
    try:
        shutil.copy(file_dir+filename+'.py',data_dir+'\\'+filename+'.py') #Copy of the script
    except IOError:
        print 'Error saving the copy of the script'
    plot2d_1.save_png()
    plot2d_2.save_png()


    # current_source.set_current(0)
    # current_source.off()
    # Vg_source.set_voltage(0)
    # Vg_source.off()
    vna.set_status('off')
    qt.mend()
