from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
from lib.math import fit
import math
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
SSB_atom2     = qt.instruments.get('SSB_atom2')

Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)

SSB_atom.set_freq_start(4)
SSB_atom.set_freq_stop(8.)
SSB_atom.set_conversion_loss(6.)
SSB_atom.set_LO_power(15)
SSB_atom.set_band_type(-1)
SSB_atom.set_IF_frequency(0.05)

# we define the board averaging
averaging = 5e3

# we define the IF power:
# Should not be more than 0 dB.
power1 = -0. # power of the first tone
#power2 = -3.775 # for rudAtt=10
power2 = -1.
# Readout Frequency [GHz]
f_cav  = 7.028
#Qubit Frequency [GHz]
f_atom = 6.2865
#f_atom = 6.2839
print 'f_atom:', f_atom
#t_pi_o2 = 20e-9


# !V 180719
# t2 = 30e-9
t_pi_o2 = 20e-9
print 't_pi_o2: ', t_pi_o2
# floor(t2*1e9/2)*1e-9 # we define the excitation time for a pi over 2 pulse on the qubit
# print 't_pi_o2=', t_pi_o2*1e9

# we define the measurement times
t_meas = 500e-9
delta_t = 200 # not to touch for now
acq_time = t_meas*1e9 + delta_t + 300. # not to touch for now

# We define the time vector of the Ramsey measurement:
tr_stop = 1e-6*4 # in s
tr_step = 5e-9*4 # in s
tr_start = 0e-6 # in s
T_vec = np.arange(tr_start, tr_stop, tr_step)


nb_sequences = len(T_vec) # we define the number of sequences to acquire for the board.

t_wait = 0e-9
t_rise = None
###########################################################
#
#
#               Experiment
#
#
###########################################################

if Tabor_loading:
    # here we write in the AWG memory if the condition Tabor_loading is True
    print 'Tabor loading...'
    Pulsing_instrument.set_trigger_time(100.)
    Pulsing_instrument.write_Ramsey_pulsessequence(t_pi_o2, tr_stop, tr_step,
            tr_start, t_meas, t_wait=t_wait, delete='all', t_rise =t_rise )


qt.mstart()

data_measurement = qt.Data(name='Ramsey')
data_measurement.add_coordinate('waiting time [ns]', units = 'ns')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Re ',            units = 'Volt')
data_measurement.add_value('Im ',            units = 'Volt')
data_measurement.create_file()


plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 ramsey',
                  coorddim  = 0,
                  valdim    = 1,
                  maxtraces = 2)

plot2d_2 = qt.Plot2D(data_measurement,
                    name      = 'Phase ramsey',
                    coorddim  = 0,
                    valdim    = 2,
                    maxtraces = 2)
#
plot2d_3 = qt.Plot2D(data_measurement,
                  name      = 'Re ramsey',
                  coorddim  = 0,
                  valdim    = 3,
                  maxtraces = 2)

plot2d_4 = qt.Plot2D(data_measurement,
                    name      = 'Im ramsey',
                    coorddim  = 0,
                    valdim    = 4,
                    maxtraces = 2)

if FIT:

    data_fit = qt.Data(name='Ramsey_OSC_fit')
    data_fit.add_value('parameters ',            units = 'rad, rad, GHz*2pi, rad, ns')
    data_fit.add_value('errors ',            units = 'rad, rad, GHz*2pi, rad, ns')
    data_fit.create_file()

try:
    Pulsing_instrument.prep_rabi(f_cav, f_atom, averaging, nb_sequences,
        power1, power2, acq_time, t_meas*1e9, delta_t)
    # Pulsing_instrument.prep_ramsey(f_cav, f_atom, averaging, nb_sequences, power1, power2)
    qt.msleep(3)
    Tabor.set_trigger_source('TIM')
    while Pulsing_instrument.get_acquisition_completed() != 100.:
        print  Pulsing_instrument.get_acquisition_completed(), '%'

        result = Pulsing_instrument.measurement()
        ((real_a, rea0), (imag_a, ima0))= result
        real_a -= rea0
        imag_a -= ima0
        amplitude = np.sqrt(real_a**2+imag_a**2)

        complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        phase = np.angle(complexe)

        # (real_a, imag_a)= result
        #
        # amplitude = np.sqrt(real_a**2+imag_a**2)
        #
        # complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        # phase = np.angle(complexe)

        qt.msleep(0.1)

        # plot2d_1.replace_inline_data(T_vec*1e9, amplitude)
        plot2d_2.replace_inline_data(T_vec*1e9, phase)
        # plot2d_3.replace_inline_data(T_vec*1e9, real_a)
        # plot2d_4.replace_inline_data(T_vec*1e9, imag_a)
        if FIT:
            s = fit.ExponentialDecaySine()
            signal = real_a
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 4.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 1000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)
            plot2d_3.replace_inline_data_y2(T_vec*1e9, signal, values_from_fit )
            plot2d_3.set_plottitle('T2= '+ str(p[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')

            s = fit.ExponentialDecaySine()
            signal = imag_a
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 10.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 1000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)
            plot2d_4.replace_inline_data_y2(T_vec*1e9, signal, values_from_fit )
            plot2d_4.set_plottitle('T2= '+ str(p[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')

            s = fit.ExponentialDecaySine()
            signal = amplitude
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 10.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            # print pulsation/2/np.pi
            phio = 0.
            decaytime = 1000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)
            # print 'params:', s.get_fit_params()
            # print 'errors:', s.get_fit_errors()

            plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )

            # plot2d_1.add(T_vec*1e9, values_from_fit)
            plot2d_1.set_plottitle('T2= '+ str(p[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')

    Tabor.set_trigger_source('EVEN')
    if Pulsing_instrument.get_board_flag():
        Pulsing_instrument.measurement_close(transfert_info=False)

finally:
    if Pulsing_instrument.get_board_flag():
        Pulsing_instrument.measurement_close(transfert_info=False)

    data_measurement.add_data_point(T_vec*1e9, amplitude, phase, real_a, imag_a)
    data_measurement.close_file()

    print Pulsing_instrument.measurement_close(transfert_info=True)
    # smb_atom.set_freqsweep('OFF')
    # smb_atom.set_gui_update('ON')
    Tabor.set_trigger_source('EVEN')

    if FIT:
        plot2d_1.add(T_vec*1e9, values_from_fit)
        data_fit.add_data_point( s.get_fit_params(), s.get_fit_errors())
        data_fit.close_file()


plot2d_1.save_png()
# plot2d_2.save_png()
qt.mend()
