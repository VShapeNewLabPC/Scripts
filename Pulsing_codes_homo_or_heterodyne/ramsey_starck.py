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


averaging = 5e3
power1 = -0 #in [dB]
power2 = -0.3
# Frequency [GHz]
f_cav  = 7.88
f_atom = 6.283

t_pi_o2 = 18e-9

amplitude_photon = 0.2 # max of amplitude is 1.

t_meas = 0.5e-6
delta_t = 200
acq_time = 1200. #in ns

tr_stop = 4e-6
tr_step = 20e-9
tr_start = 0e-6

T_vec = np.arange(tr_start, tr_stop, tr_step)
nb_sequences = len(T_vec)

###########################################################
#
#
#               Experiment
#
#
###########################################################

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(100.)
    Pulsing_instrument.write_Ramsey_Starckshift_pulsessequence(t_pi_o2, tr_stop, tr_step, tr_start, amplitude_photon, t_meas, t_protect=0, delete='all')
    # print 'here'

qt.mstart()

data_measurement = qt.Data(name='Ramsey')
data_measurement.add_coordinate('waiting time [ns]', units = 'ns')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Re ',            units = 'Volt')
data_measurement.add_value('Im ',            units = 'Volt')
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

if FIT:
    data_fit = qt.Data(name='Ramsey_stark_OSC_fit')
    data_fit.add_value('parameters ',            units = 'rad, rad, GHz*2pi, rad, ns')
    data_fit.add_value('errors ',            units = 'rad, rad, GHz*2pi, rad, ns')
    data_fit.create_file()

board_flag = None
try:

    Pulsing_instrument.prep_rabi(f_cav, f_atom, averaging, nb_sequences,
        power1, power2, acq_time, t_meas*1e9, delta_t)  ###Same procedure for rabi and ramsey



    qt.msleep(2)

    board_flag = True

    Tabor.set_trigger_source('TIM')
    while ats9360.get_completed_acquisition() != 100.:
        print ats9360.get_completed_acquisition(), '%'
        result = ats9360.measurement()

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

        if FIT:
            s = fit.ExponentialDecaySine()
            signal = amplitude
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 10.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 2000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            p_err = s.get_fit_errors()
            values_from_fit = s.func(p)
            # print 'params:', s.get_fit_params()
            # print 'errors:', s.get_fit_errors()
            plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )
            plot2d_1.set_plottitle('T2= '+ str(p[4])+' +- '+str(p_err[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' +- ' +str(p_err[2]/2/np.pi*1e3) +' MHz')

            signal = phase
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 10.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 2000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p_phase = s.fit(p0)
            values_from_fit_phase = s.func(p_phase)
            p_phase_err = s.get_fit_errors()
            plot2d_2.replace_inline_data_y2(T_vec*1e9, phase, values_from_fit_phase)
            plot2d_2.set_plottitle('T2= '+ str(p_phase[4])+' +- '+str(p_phase_err[4])+' ns'+', Freq= '+str(p_phase[2]/2/np.pi*1e3)+' MHz')
            signal = real_a
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 10.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 2000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            p_err = s.get_fit_errors()
            values_from_fit = s.func(p)
            plot2d_3.replace_inline_data_y2(T_vec*1e9, real_a, values_from_fit )
            plot2d_3.set_plottitle('T2= '+ str(p[4])+' +- '+str(p_err[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')
            signal = imag_a
            s.set_data(T_vec*1e9, signal)
            # guess parameters##########################################################
            background = (signal.max() + signal.min() )/2.
            osc_amp = (signal.max() - signal.min() )/2.
            nb_expected_oscillation = 10.
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 2000.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            p_err = s.get_fit_errors()
            values_from_fit = s.func(p)
            plot2d_4.replace_inline_data_y2(T_vec*1e9, imag_a, values_from_fit )
            plot2d_4.set_plottitle('T2= '+ str(p[4])+' +- '+str(p_err[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')
            # plot2d_4.replace_inline_data(T_vec*1e9, imag_a)

        else:
            plot2d_1.replace_inline_data(T_vec*1e9, amplitude)
            plot2d_2.replace_inline_data(T_vec*1e9, phase)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)

    data_measurement.add_data_point(T_vec*1e9,amplitude,phase, real_a, imag_a)


    print ats9360.measurement_close(transfert_info=True)


    Tabor.set_trigger_source('EVEN')

plot2d_1.save_png()
plot2d_2.save_png()
if FIT:

    data_fit.add_data_point( p,p_err)
    data_fit.add_data_point( p_phase,p_phase_err)
    data_fit.close_file()
    plot2d_1.add(T_vec*1e9, values_from_fit)
    plot2d_2.add(T_vec*1e9, values_from_fit_phase)

data_measurement.close_file()


qt.mend()
