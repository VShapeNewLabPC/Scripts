from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
###########################################################
Tabor_loading = 1
###########################################################
#
#
#               Devices
#
#
###########################################################

Tabor           = qt.instruments.get('Tabor')
smb_cavity     = qt.instruments.get('smb_1')
smb_atom        = qt.instruments.get('smb_2')
ats9360        = qt.instruments.get('ats9360')
SSB_cavity      = qt.instruments.get('SSB_cavity')
SSB_atom      = qt.instruments.get('SSB_atom')
RUDAT         = qt.instruments.get('RUDAT_ph')
Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)
###########################################################
SSB_atom.set_freq_start(4)
SSB_atom.set_freq_stop(8)
SSB_atom.set_conversion_loss(6.)
SSB_atom.set_LO_power(15)
SSB_atom.set_band_type(-1)
SSB_atom.set_IF_frequency(0.05)
###########################################################

power1 = -0. #in [dB]
power2 = -0
# Frequency [GHz]
# f_pi  = 6.26922
# f_pi = 7.601
f_pi = 6.29622

#StartFrequency [GHz]
f_min = 6900e-3
# f_min = 7.4
#
#StopFrequency [GHz]
# f_max = 7.05
# f_max = 7.5
f_max = 7100e-3

#Step of the Sweep [GHz]
f_step = 0.5e-3

freq_vec = np.arange(f_min, f_max + f_step, f_step)
if len(freq_vec) %2 !=0:
        freq_vec = np.arange(f_min, f_max , f_step)

t_pi = 50e-9
# t_pi = 38e-9
#
t_ro = 0.5e-6
delta_t = 0.2e-6
acq_time = (t_ro + delta_t )*1e9 + 300
# t_rise= 10e-9
t_rise= 0.0e-9

Att_vec = np.arange(0.0, 30.0, 1.0)

av = 1e3
average_vec = np.ones_like(Att_vec)
for i in np.arange(len(Att_vec)):
    if i%2 == 0:
        # av += 100
        # av = av + 100
        av = av*2
    average_vec[i] = av
# average_vec = 100 + 100*np.round(Att_vec, 0)
t_rise = None
###########################################################
#
#
#               Experiment
#
#
###########################################################
if Tabor_loading:
    # print 'here'
    Pulsing_instrument.set_trigger_time(100.)
    Pulsing_instrument.write_twotone_pulsessequence_withpi(temp_1=t_ro,
        t1_start= t_pi + 0.1e-6, temp_2=t_pi , m1_start= t_pi ,  delete = 0, t_rise=t_rise)


# Pulsing_instrument.write_twotone_pulsessequence( 500e-9, 100e-9 + t_pi, t_pi, delete = 'all')
qt.mstart()

data_measurement = qt.Data(name='Cavity_shift')
data_measurement.add_coordinate('R.O. frequency [GHz]', units = 'GHz')
data_measurement.add_coordinate('power [dBm]', units = 'dBm')

data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Re ',            units = 'Volt')
data_measurement.add_value('Im ',            units = 'Volt')
data_measurement.add_value('S21 pi',            units = 'Volt')
data_measurement.add_value('Phase pi',            units = 'rad')
data_measurement.add_value('Re pi',            units = 'Volt')
data_measurement.add_value('Im pi',            units = 'Volt')
data_measurement.add_value('D blobs',            units = 'Volt')

data_measurement.create_file()

plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 ',
                  coorddim  = 0,
                  valdim    = 2,
                  maxtraces = 2)

# plot3d_1 = qt.Plot3D(data_measurement,
#                   name      = 'S21 colorplot ',
#                   coorddim  = (0,1),
#                   valdim    = 2)

plot2d_pi = qt.Plot2D(data_measurement,
                  name      = 'S21pi ',
                  coorddim  = 0,
                  valdim    = 6,
                  maxtraces = 2)

plot2d_d = qt.Plot2D(data_measurement,
                  name      = 'd blobs ',
                  coorddim  = 0,
                  valdim    = 10,
                  maxtraces = 2)
# plot3d_pi = qt.Plot3D(data_measurement,
#                   name      = 'S21pi colorplot ',
#                   coorddim  = (0,1),
#                   valdim    = 6)
# plot2d_2 = qt.Plot2D(data_measurement,
#                     name      = 'Phase ',
#                     coorddim  = 0,
#                     valdim    = 2,
#                     maxtraces = 2)
# #
# plot2d_3 = qt.Plot2D(data_measurement,
#                   name      = 'Re ',
#                   coorddim  = 0,
#                   valdim    = 3,
#                   maxtraces = 2)
#
# plot2d_4 = qt.Plot2D(data_measurement,
#                     name      = 'Im ',
#                     coorddim  = 0,
#                     valdim    = 4,
#                     maxtraces = 2)
plot3d_1 = qt.Plot3D(data_measurement,
                  name      = 'd blobs colorplot ',
                  coorddim  = (0,1),
                  valdim    = 10)


for i, att in enumerate(Att_vec):
    RUDAT.set_attenuation(att)
    print att, average_vec[i]
    qt.msleep(0.1)
    board_flag = None
    try:
        Pulsing_instrument.prep_conditional_transmission(freq_vec, average_vec[i],
                    power1, f_cw=f_pi, power2=power2, acq_time=acq_time, pulse_time=t_ro*1e9, delta_t=delta_t )
        # prep_shift_homodyne(acq_time, delta_t, freq_vec, averaging, power1, 'two' , f_pi, power2 )
        # qt.msleep(1)
        smb_atom.set_freqsweep('OFF')
        smb_cavity.set_freqsweep('ON')
        smb_cavity.restartsweep()

        qt.msleep(2)

        board_flag = True
        Tabor.set_trigger_source('TIM')
        while ats9360.get_completed_acquisition() != 100.:
            print  ats9360.get_completed_acquisition(), '%'
            result = ats9360.measurement()
            # (real, imag)= result
            ((real, rea0), (imag,ima0))= result
            real = real - rea0
            imag = imag - ima0

            real = np.reshape(real, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) )
            imag = np.reshape(imag, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) )
            real_a = real[:, :Pulsing_instrument.get_pulsenumber_averaging()/2]
            imag_a = imag[:, :Pulsing_instrument.get_pulsenumber_averaging()/2]
            real_a = np.mean(real_a, axis = 1)
            imag_a = np.mean(imag_a, axis = 1)
            real_a_pi = real[:, Pulsing_instrument.get_pulsenumber_averaging()/2:]
            imag_a_pi = imag[:, Pulsing_instrument.get_pulsenumber_averaging()/2:]
            real_a_pi = np.mean(real_a_pi, axis = 1)
            imag_a_pi = np.mean(imag_a_pi, axis = 1)

            d_blobs = np.sqrt((real_a-real_a_pi)**2+(imag_a-imag_a_pi)**2)

            amplitude = np.sqrt(real_a**2+imag_a**2)
            complexe = (real_a + 1j*imag_a )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase = np.angle(complexe)

            amplitude_pi = np.sqrt(real_a_pi**2+imag_a_pi**2)
            complexe_pi = (real_a_pi + 1j*imag_a_pi )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase_pi = np.angle(complexe_pi)
            qt.msleep(0.1)
            #
            plot2d_1.replace_inline_data(freq_vec, amplitude)
            plot2d_pi.replace_inline_data(freq_vec, amplitude_pi)
            plot2d_d.replace_inline_data(freq_vec, d_blobs)


        Tabor.set_trigger_source('EVEN')
        ats9360.measurement_close(transfert_info=False)
        board_flag = False
    finally:
        if board_flag:
            ats9360.measurement_close(transfert_info=False)

        data_measurement.add_data_point(freq_vec,att*np.ones_like(freq_vec),
            amplitude,phase, real_a, imag_a,
            amplitude_pi, phase_pi, real_a_pi, imag_a_pi,
            d_blobs)
        data_measurement.new_block()
        print ats9360.measurement_close(transfert_info=True)
        # smb_cavity.set_freqsweep('OFF')
        # smb_cavity.set_gui_update('ON')
        Tabor.set_trigger_source('EVEN')


data_measurement.close_file()
plot3d_1.save_png()
# plot3d_pi.save_png()
# plot2d_3.save_png()
# plot2d_4.save_png()
qt.mend()
