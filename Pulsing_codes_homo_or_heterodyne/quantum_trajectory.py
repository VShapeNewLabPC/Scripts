from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
import matplotlib.pyplot as plt

Tabor_loading=1
###########################################################
#
#
#               Devices
#
#
###########################################################
Tabor           = qt.instruments.get('Tabor')
smb_cavity     = qt.instruments.get('smb_1')
ats9360        = qt.instruments.get('ats9360')
SSB_cavity      = qt.instruments.get('SSB_cavity')
Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.)

averaging = 500
power1 = -0.0
power2 = -1
t2 = 40e-9 # should be a pi pulse
freq_cav = 7.028
##(f_atom=f_qubit)
freq_ex_pi = 6.2865

t1 = 600e-9 #duration of readout pulse
t1_1 = 50e-9 #?for first readout in measurment (g)

t1_2 = t1 #?for second readout in measurment (g or e)
# 30.
delta_t = 200.

# acq_time = 2600.
t1_1start = 0.1e-6
t_between = 300e-9
t2_start = t1_1start + t1_1 + t_between
t1_2start = t2_start + t2

t_protect = 0.
delay_initial = 100.
ats9360.set_trigger_delay(delay_initial)
acq_time = (t1_1start+t1)*1e9+300

delta_m1_start = 0.1e-6
###########################################################
#
#
#               Experiment
#
#
###########################################################

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(200.)
    Pulsing_instrument.write_IQ_alwayspi_several_RO(0., t1_1, t1_2, t2_start=t2_start ,
            t1_1start=t1_1start, t1_2start=t1_2start, delta_m1_start=0.1e-6, delete=False)


ats9360.set_acquisition_time(acq_time)


qt.mstart()
#
data_measurement = qt.Data(name='Quantum_time')
data_measurement.add_coordinate('time ',         units = 'us')
data_measurement.add_value('I (t) ',         units = 'Volt')
data_measurement.add_value('Q (t) ',         units = 'Volt')
data_measurement.add_value('Ipi (t) ',         units = 'Volt')
data_measurement.add_value('Q pi(t) ',         units = 'Volt')
data_measurement.create_file()

# without pi pulse

board_flag = None
try:
    Pulsing_instrument.prep_timing( freq_cav, averaging, power1, average_type='test')
    qt.msleep(2)
    # commented by remy 20181121
    # Tabor.set_ch4_output('ON')
    # amplitude2 = 10**((power2)/10.)
    # Tabor.set_ch4_amplitude(2*amplitude2)
    # print 'here', Tabor.get_ch4_output(), Tabor.get_ch4_amplitude()
    board_flag = True
    Tabor.set_trigger_source('TIM')
    while ats9360.get_completed_acquisition() != 100.:
        print 'percentage done:' + str(ats9360.get_completed_acquisition())+'%'

        result = ats9360.measurement()
        (I, Q) = result
        qt.msleep(0.1)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)
    print ats9360.measurement_close(transfert_info=True)
Tabor.set_trigger_source('EVEN')

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(200.)
    Pulsing_instrument.write_IQ_alwayspi_several_RO(t2, t1_1, t1_2, t2_start=t2_start ,
            t1_1start=t1_1start, t1_2start=t1_2start, delta_m1_start=0.1e-6, delete=False)
# with pi pulse
ats9360.set_trigger_delay(delay_initial + (t1_1+t_between+t2)*1e9 )
board_flag = None
try:
    Pulsing_instrument.prep_timing( freq_cav, averaging, power1, average_type='test')
    qt.msleep(2)
    Tabor.set_ch4_output('ON')
    amplitude2 = 10**((power2)/10.)
    Tabor.set_ch4_amplitude(2*amplitude2)
    # print 'here', Tabor.get_ch4_output(), Tabor.get_ch4_amplitude()
    board_flag = True
    Tabor.set_trigger_source('TIM')
    while ats9360.get_completed_acquisition() != 100.:
        print 'percentage done:' + str(ats9360.get_completed_acquisition())+'%'

        result = ats9360.measurement()
        (Ipi, Qpi) = result
        qt.msleep(0.1)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)
    print ats9360.measurement_close(transfert_info=True)
Tabor.set_trigger_source('EVEN')


for i in np.arange(len(I[:,0])):
    data_measurement.add_data_point(np.arange(len(I[0,:]))\
                /ats9360.get_samplerate()*1e3, I[i,:], Q[i,:], Ipi[i,:], Qpi[i,:])
    data_measurement.new_block()


data_measurement.close_file()

data_stat = qt.Data(name='time_plot')
data_stat.add_coordinate('time ',         units = 'us')
data_stat.add_value('Imean ',         units = 'Volt')
data_stat.add_value('Qmean ',         units = 'Volt')
data_stat.add_value('Ipimean ',         units = 'Volt')
data_stat.add_value('Qpimean ',         units = 'Volt')
data_stat.add_value('Istd ',         units = 'Volt')
data_stat.add_value('Qstd ',         units = 'Volt')
data_stat.add_value('Ipistd ',         units = 'Volt')
data_stat.add_value('Qpistd ',         units = 'Volt')
data_stat.create_file()

Imean = np.zeros(len(I[0,:]))
Qmean = np.zeros_like(Imean)
Ipimean = np.zeros_like(Imean)
Qpimean = np.zeros_like(Imean)
Istd = np.zeros_like(Imean)
Qstd = np.zeros_like(Imean)
Ipistd = np.zeros_like(Imean)
Qpistd = np.zeros_like(Imean)
time_vec = np.arange(len(I[0,:]))/ats9360.get_samplerate()*1e3
for i in np.arange(len(Imean)):
    Imean[i] = np.mean(I[:,i])
    Ipimean[i] = np.mean(Ipi[:,i])
    Qmean[i] = np.mean(Q[:,i])
    Qpimean[i] = np.mean(Qpi[:,i])
    Istd[i] = np.std(I[:,i])
    Ipistd[i] = np.std(Ipi[:,i])
    Qstd[i] = np.std(Q[:,i])
    Qpistd[i] = np.std(Qpi[:,i])

data_stat.add_data_point(np.arange(len(Imean)), Imean, Qmean, Ipimean, Qpimean,
                                                 Istd, Qstd, Ipistd, Qpistd)

data_stat.close_file()
qt.mend()

fig, ax = plt.subplots(1,1)
ax.grid()
ax.plot(time_vec, Imean, '-', color='b', linewidth=1.5)
ax.plot(time_vec, Ipimean, '-', color='r', linewidth=1.5)

ax.plot(time_vec, I[0,:], '.-', color='b', alpha=0.5)
ax.plot(time_vec, Ipi[0,:], '.-', color='r', alpha=0.5)

ax.fill_between(time_vec, Imean-Istd, Imean+Istd, color='b', alpha=0.4)
ax.fill_between(time_vec, Ipimean-Ipistd, Ipimean+Ipistd, color='r', alpha=0.4)

figQ, axQ = plt.subplots(1,1)
axQ.grid()
axQ.plot(time_vec, Qmean, '-', color='b', linewidth=1.5)
axQ.plot(time_vec, Qpimean, '-', color='r', linewidth=1.5)

axQ.plot(time_vec, Q[0,:], '.-', color='b', alpha=0.5)
axQ.plot(time_vec, Qpi[0,:], '.-', color='r', alpha=0.5)

axQ.fill_between(time_vec, Qmean-Qstd, Qmean+Qstd, color='b', alpha=0.4)
axQ.fill_between(time_vec, Qpimean-Qpistd, Qpimean+Qpistd, color='r', alpha=0.4)
plt.show()
