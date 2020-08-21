from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
from lib.math import fit
###########################################################
Tabor_loading = 0
Tabor_loading = 1

FIT = False
FIT = True
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

Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)

SSB_atom.set_freq_start(4)
SSB_atom.set_freq_stop(8)
SSB_atom.set_conversion_loss(6.)
SSB_atom.set_LO_power(15)
SSB_atom.set_band_type(-1)
SSB_atom.set_IF_frequency(0.05)

averaging = 2e3
power1 = -0. #in [dB]
power2 = -0.
# Frequency [GHz]
f_cav  = 7.535

#Atom Frequency [GHz]
f_atom = 6.2612

# t_pi = 110e-9
t_pi = 36e-9


t_meas = 500e-9
delta_t = 200.
acq_time = 1200.

T_vec = np.logspace(-9., -4.3, 100)

before = 10
t_prev = np.linspace(-T_vec[-1]/5, 0, before)

amplitude_RO = 1.
amplitude_photon = 1
# print amplitude_photon
###########################################################
#
#
#               Experiment
#
#
###########################################################

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(500.)
    Pulsing_instrument.write_Relaxation_pulsessequence_with_photons(t_pi, T_vec, amplitude_photon, amplitude_RO,  t_meas, delete='all', before=before)


T_vec = np.append(t_prev, T_vec)
nb_sequences = len(T_vec)
qt.mstart()

data_measurement = qt.Data(name='Relaxation_n=' + str(amplitude_photon))
data_measurement.add_coordinate('excitation time [ns]', units = 'ns')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Re ',            units = 'Volt')
data_measurement.add_value('Im ',            units = 'Volt')
data_measurement.add_value(' log',            units = '')

data_measurement.create_file()


plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 ',
                  coorddim  = 0,
                  valdim    = 1,
                  maxtraces = 2)

plot2d_2 = qt.Plot2D(data_measurement,
                    name      = 'Phase ',
                    coorddim  = 0,
                    valdim    = 2,
                    maxtraces = 2)

plot2d_3 = qt.Plot2D(data_measurement,
                  name      = 'Re ',
                  coorddim  = 0,
                  valdim    = 3,
                  maxtraces = 2)

plot2d_4 = qt.Plot2D(data_measurement,
                    name      = 'Im ',
                    coorddim  = 0,
                    valdim    = 4,
                    maxtraces = 2)

plot2d_5  = qt.Plot2D(data_measurement,
                    name      = 'log amp ',
                    coorddim  = 0,
                    valdim    = 5,
                    maxtraces = 2)
if FIT:
    data_fit = qt.Data(name='T1_fit')
    data_fit.add_value('parameters ',            units = 'rad, rad, -GHz, .. ')
    data_fit.add_value('perrors ',            units = 'rad, rad, -GHz, ..')
    data_fit.create_file()
board_flag = None
try:
    Pulsing_instrument.prep_relaxation(f_cav, f_atom, averaging, nb_sequences,
        power1, power2, acq_time, t_meas*1e9, delta_t)
    qt.msleep(2.)

    board_flag = True

    Tabor.set_trigger_source('TIM')
    while ats9360.get_completed_acquisition() != 100.:
        print ats9360.get_completed_acquisition()
        result = ats9360.measurement()
        ((real_a, rea0), (imag_a, ima0))= result
        real_a -= np.mean(rea0)
        imag_a -= np.mean(ima0)
        amplitude = np.sqrt(real_a**2+imag_a**2)

        complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        phase = np.angle(complexe)
        # (real_a, imag_a)= result

        # real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        # imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        # amplitude = np.sqrt(real_a**2+imag_a**2)
        #
        # complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        # phase = np.angle(complexe)
        # phase_renormed = np.log( phase[-10:-1].mean() - phase)
        # phase_renormed = np.log(phase_renormed)
        signal = imag_a
        a0 = np.mean(signal[:before])
        log_amp = np.log(a0 - signal)
        qt.msleep(0.1)

        if FIT:
            s = fit.Exponential()

            a0 = np.mean(signal[:before])
            sig = a0 - signal

            s.set_data(T_vec[before:]*1e9, sig[before:] )
            # guess parameters##########################################################
            a = 0.
            b = np.max(sig)
            c = 0.
            d = 1./2e3
            p0 = [a, b, c, d]

            # fitting ##################################################################
            p = s.fit(p0, fixed = [0, 2])
            # p = s.fit(p0)
            values_from_fit = s.func(p)
            # print 'params:', s.get_fit_params()
            # print 'errors:', s.get_fit_errors()
            # print 'T1:', 1/p[3], 'ns'
            values_from_fit= np.append( a0+np.zeros(before),  a0 - values_from_fit)
            plot2d_4.replace_inline_data_y2(T_vec*1e9, signal, values_from_fit)
            plot2d_4.set_plottitle('T1= '+ str(1/p[3])+' ns')

        else :
            plot2d_4.replace_inline_data(T_vec*1e9, signal)

        plot2d_1.replace_inline_data(T_vec*1e9, amplitude)
        plot2d_2.replace_inline_data(T_vec*1e9, phase)
        plot2d_3.replace_inline_data(T_vec*1e9, real_a)
        # plot2d_4.replace_inline_data(T_vec*1e9, imag_a)
        # plot2d_5.replace_inline_data(T_vec*1e9, log_amp)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)

    data_measurement.add_data_point(T_vec*1e9,amplitude,phase, real_a, imag_a, log_amp)
    data_measurement.close_file()

    print ats9360.measurement_close(transfert_info=True)
    # smb_atom.set_freqsweep('OFF')
    # smb_atom.set_gui_update('ON')
    Tabor.set_trigger_source('EVEN')

if FIT:
    data_fit.add_data_point(s.get_fit_params(), s.get_fit_errors() )
    data_fit.close_file()
    plot2d_4.add(T_vec*1e9, values_from_fit)

    # plot2d_5.add(T_vec*1e9, np.log(np.abs(values_from_fit)) )
    plot2d_4.set_plottitle('T1= '+ str(1/p[3])+' ns')
    # plot2d_5.set_plottitle('T1= '+ str(1/p[3])+' ns')



plot2d_1.save_png()
plot2d_2.save_png()
plot2d_4.save_png()

# plot2d_5.save_png()

qt.mend()
