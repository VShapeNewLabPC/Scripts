# -*- coding: utf-8 -*-
import qt
import numpy as np
import shutil

filename = 'VNA_vs_B_standard' #Filename that will be used to name the data folder. Use script filename.
file_dir = '' #File directory. It has the form working_directory\file_dir\filename

#We create the object data
data = qt.Data(name=str(filename)) #Do not change filename. Include comments as 'filename + str(comments)'
data_dir=qt.Data.get_dir(data)           #Where to save the copy of the script.

try:
    shutil.copy(file_dir+filename+'.py',data_dir+'\\'+filename+'.py') #Copy of the script
except IOError:
    print 'Error saving the copy of the script'
###########################################################
#
#
# Sweep and measurement parameters
#
#
###########################################################

#Number of points (int) [Sample]
points = 101

#StartFrequency [Hz]
startfrequency =6e9

#StopFrequency [Hz]
stopfrequency =8e9

#Power [dBm]
detector_power = -30.# dBm

#Detector bandwidth [Hz]
detector_bandwidth = 10000 # Hz

#Averages [Sample]
averages = 1

#Sweeps
sweeps = averages



#Coil current sweep: start, stop, step, current for reference measurement
c_start = -3e-3 #Ampere
c_stop  = 3e-3
c_step  = 1e-3
c_ref = 0e-3
curr_vec = np.arange(c_start,c_stop + c_step, c_step)
#curr_vec = np.array([20.19e-3,18.48e-3,18.61e-3,18.79e-3,19.181e-3])
#Sparameter
Sparam='S21'


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


znb = qt.instruments.get('znb')

current_source = qt.instruments.get('hp3245')
#Vg_source = qt.instruments.get('hp3245a_2')

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


znb.initialize_one_tone_spectroscopy(('trace1',),(Sparam,))


znb.set_startfrequency(startfrequency)
znb.set_stopfrequency(stopfrequency)
znb.set_averages(averages)
znb.set_averagestatus('on')
znb.set_sweeps(sweeps)
znb.set_power(detector_power)
znb.set_measBW(detector_bandwidth)
znb.set_points(points)

###########################################################
#
#
#                      HP3245A
#
#
###########################################################

current_source.set_output_terminal('FRONT') #We use the rear ports
current_source.set_mode('dci')
current_source.set_channel('A')
current_source.set_resolution('high')


###########################################################
#
#
#               Experiment
#
#
###########################################################

qt.mstart()



#We add the coordinate of what is swept
data.add_coordinate('Current [mA]', units='mA')
data.add_coordinate('Frequency [GHz]', units='GHz')

#We add the value that will be read out

data.add_value('S_21 [dB]', units='dB')
data.add_value('Phase [deg]', units='deg')
data.add_value('S_21_normalized [dB]', units='dB')
data.add_value('Phase_normalized [deg]', units='deg')

#We live plot what we record

plot2d_1= qt.Plot2D(data, name='last_trace', coorddim=1, valdim=2, maxtraces=2)
plot2d_2= qt.Plot2D(data, name='last_trace_phase', coorddim=1, valdim=3, maxtraces=2)


plot3d_1 = qt.Plot3D(data, name='S21_dB', coorddims=(0,1), valdim=2,style='image')
plot3d_2 = qt.Plot3D(data, name='S21_phase', coorddims=(0,1), valdim=3,style='image')
plot3d_3 = qt.Plot3D(data, name='S21_dB_normalized', coorddims=(0,1), valdim=4,style='image')
plot3d_4 = qt.Plot3D(data, name='S21_phase_normalized', coorddims=(0,1), valdim=5,style='image')
plot3d_1.set_palette('bluewhitered')
plot3d_2.set_palette('bluewhitered')
plot3d_3.set_palette('bluewhitered')
plot3d_4.set_palette('bluewhitered')

#We create a file corresponding to the object data
data.create_file()

# Vg_source.on()
# Vg_source.set_voltage(Vg_ref)

current_source.on()
qt.msleep(2)
current_source.set_current(c_ref)
qt.msleep(2)

znb.averageclear()
# znb.measure()
qt.msleep(0.1)
#ref_amp,ref_phase=znb.get_traces( ('trace1',))[0]
ref_amp = 0
ref_phase = 0
qt.msleep(0.1)
znb.set_averages(averages)
freqs=linspace(startfrequency,stopfrequency,points)/1.0e9

try:
    for c in curr_vec:
		current_source.set_current(c)
		qt.msleep(1)
		znb.averageclear()
		znb.measure()
		qt.msleep(0.1)
		amp,phase=znb.get_traces( ('trace1',))[0]
		qt.msleep(0.1)
		amp_norm = amp - ref_amp
		phase_norm= phase - ref_phase
        #We save data
		data.add_data_point(c*np.ones_like(freqs)*1e3,freqs,amp,phase,amp_norm,phase_norm)
		data.new_block()
except Exception as error:
    print str(error)
finally:
    qt.msleep()
    data.close_file()
    plot3d_1.save_png()
    plot3d_2.save_png()
    plot3d_3.save_png()
    plot3d_4.save_png()

    current_source.set_current(0)
    current_source.off()
    # Vg_source.set_voltage(0)
    # Vg_source.off()
    znb.set_status('off')
    qt.mend()
