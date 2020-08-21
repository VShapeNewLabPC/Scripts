from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
import datetime
###########################################################
Tabor_loading = 0
# Tabor_loading = 1

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
# SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)

SSB_atom.set_freq_start(4)
SSB_atom.set_freq_stop(8.)
SSB_atom.set_conversion_loss(6.)
SSB_atom.set_LO_power(15)
SSB_atom.set_band_type(-1)
SSB_atom.set_IF_frequency(0.05)

# SSB_3.set_freq_start(4)
# SSB_3.set_freq_stop(8.)
# SSB_3.set_conversion_loss(6.)
# SSB_3.set_LO_power(15)
# SSB_3.set_band_type(-1)
# SSB_3.set_IF_frequency(0.05)
###########################################################

averaging = 2e3
#amplitude of modulated signal (==voltage on mixers)
power1 = -0.0 #in [dB]
power2 = -1.
# -0.55

# Frequency [GHz]
f_cav = 7.028
#f_cav = 6.95
#Atom Frequency [GHz]
f_atom = 6.2865
# time of readout pulse
t_meas = 500e-9

delta_t = 200
acq_time = t_meas*1e9 + delta_t + 300.

tr_stop = 0.2e-6 #increased 2 times - ok
tr_step = 1e-9
#tr_stop = 2e-6/10
#tr_step = 10e-9/10
nom_expect = 1.

tr_start = 0e-9
t_wait = 0e-9

T_vec = np.arange(tr_start, tr_stop, tr_step)
nb_sequences = len(T_vec)

t_rise = None
# t_rise = 10e-9
###########################################################
#               Save Parameters
print 'pwr1:', power1, ' pwr2:', power2, ' f_cav:', f_cav,' f_atom:', f_atom,' t_meas:', t_meas
qub_cur_sour = qt.instruments.get('hp3245')
print 'current =', qub_cur_sour.get_current()

RUDAT = qt.instruments.get('RUDAT_ph')
print 'rudAtt', RUDAT.get_attenuation()

###########################################################
#
#
#               Experiment
#
#
###########################################################
mw = 2
# Pulsing_instrument.set_routing_awg({'secondtone_channel':4})

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(100.)
    Pulsing_instrument.write_Rabi_pulsessequence(tr_stop, tr_step, tr_start, t_meas,
            t_wait=t_wait, delta_m1_start=0.1e-6, phi=0, delete='all', t_rise =t_rise )
    # Pulsing_instrument.write_RabiGaussian_pulsessequence(tr_stop, tr_step, tr_start, t_meas,
    #         t_wait=t_wait, delta_m1_start=0.1e-6, phi=0, delete='all', t_rise =t_rise, Nsigma=3. )


qt.mstart()

data_measurement = qt.Data(name='Rabi')
data_measurement.add_coordinate('excitation time [ns]', units = 'ns')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Re ',            units = 'Volt')
data_measurement.add_value('Im ',            units = 'Volt')
data_measurement.create_file()


plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 rabi',
                  coorddim  = 0,
                  valdim    = 1)

plot2d_2 = qt.Plot2D(data_measurement,
                    name      = 'Phase rabi',
                    coorddim  = 0,
                    valdim    = 2,
                    maxtraces = 2)
#
plot2d_3 = qt.Plot2D(data_measurement,
                  name      = 'Re rabi',
                  coorddim  = 0,
                  valdim    = 3,
                  maxtraces = 2)

plot2d_4 = qt.Plot2D(data_measurement,
                    name      = 'Im rabi',
                    coorddim  = 0,
                    valdim    = 4,
                    maxtraces = 2)
if FIT:
    data_fit = qt.Data(name='Rabi_OSC_fit')
    data_fit.add_value('parameters ',            units = 'rad, rad, GHz*2pi, rad, ns')
    data_fit.add_value('errors ',            units = 'rad, rad, GHz*2pi, rad, ns')
    data_fit.create_file()
    # TT_vec = np.append(T_vec, T_vec)
# board_flag = None
try:
    Pulsing_instrument.prep_rabi(f_cav, f_atom, averaging, nb_sequences,
        power1, power2, acq_time, t_meas*1e9, delta_t, mw =mw)
    qt.msleep(2.1)

    Tabor.set_trigger_source('TIM')
    while Pulsing_instrument.get_acquisition_completed() !=100.:
        print  Pulsing_instrument.get_acquisition_completed(), '%'

        result = Pulsing_instrument.measurement()
        ((real_a, rea0), (imag_a, ima0))= result
        real_a -= rea0
        imag_a -= ima0
        amplitude = np.sqrt(real_a**2+imag_a**2)

        complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        phase = np.angle(complexe)

        qt.msleep(0.1)

        # plot2d_1.replace_inline_data(T_vec*1e9, amplitude)
        plot2d_2.replace_inline_data(T_vec*1e9, phase)
        # plot2d_3.replace_inline_data(T_vec*1e9, real_a)
        # plot2d_4.replace_inline_data(T_vec*1e9, imag_a)
        if FIT:
            # s = fit.Sine()
            # s.set_data(T_vec*1e9, amplitude)
            # # guess parameters##########################################################
            # background = (amplitude.max() + amplitude.min() )/2.
            # osc_amp = (amplitude.max() - amplitude.min() )/2.
            # nb_expected_oscillation = 5
            # pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            # phio = 0.
            # p0 = [background, osc_amp, pulsation, phio]
            # # fitting ##################################################################
            # p = s.fit(p0)
            # values_from_fit = s.func(p)
            # # print 'params:', s.get_fit_params()
            # # print 'errors:', s.get_fit_errors()
            # plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )
            # # plot2d_1.add(T_vec*1e9, values_from_fit)
            # plot2d_1.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns')

            s = fit.ExponentialDecaySine()
            s.set_data(T_vec*1e9, amplitude)
            # guess parameters##########################################################
            background = (amplitude.max() + amplitude.min() )/2.
            osc_amp = (amplitude.max() - amplitude.min() )/2.
            nb_expected_oscillation = nom_expect
            pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
            phio = 0.
            decaytime = 1000.
            # print pulsation/2/np.pi

            #----!V added by V 180823-----
            #_Looking for number of oscillation in time
            Find_nom_expect = True
            if Find_nom_expect:
                pulsation_arr = list()
                errors_arr = list()
                pulsation_arr.append(0)
                for i in range(1,20):
                    pulsation_arr.append(2*np.pi/(tr_stop*1e9)*i)
                    p0 = [background, osc_amp, pulsation_arr[i], phio, decaytime]
                    p = s.fit(p0)
                    values_from_fit = s.func(p)
                    dif = amplitude - values_from_fit
                    errors_arr.append(sum(abs(dif)))

                id = argmin(errors_arr)
                pulsation = pulsation_arr[id]
            #----!V_end--------------------



            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)
            values = values_from_fit

            # print 'params:', s.get_fit_params()
            # print 'errors:', s.get_fit_errors()
            plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )
            # plot2d_1.add(T_vec*1e9, values_from_fit)
            plot2d_1.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')
            # print 'pi pulse time is : '+ str(np.pi/p[2]) +' ns'

            s = fit.ExponentialDecaySine()
            s.set_data(T_vec*1e9, imag_a)
            # guess parameters##########################################################
            background = (imag_a.max() + imag_a.min() )/2.
            osc_amp = (imag_a.max() - imag_a.min() )/2.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)

            plot2d_4.replace_inline_data_y2(T_vec*1e9, imag_a, values_from_fit )
            plot2d_4.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')

            s = fit.ExponentialDecaySine()
            s.set_data(T_vec*1e9, real_a)
            # guess parameters##########################################################
            background = (real_a.max() + real_a.min() )/2.
            osc_amp = (real_a.max() - real_a.min() )/2.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)

            plot2d_3.replace_inline_data_y2(T_vec*1e9, real_a, values_from_fit )
            plot2d_3.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')

            s = fit.ExponentialDecaySine()
            s.set_data(T_vec*1e9, phase)
            # guess parameters##########################################################
            background = (phase.max() + phase.min() )/2.
            osc_amp = (phase.max() - phase.min() )/2.
            p0 = [background, osc_amp, pulsation, phio, decaytime]
            # fitting ##################################################################
            p = s.fit(p0)
            values_from_fit = s.func(p)

            plot2d_2.replace_inline_data_y2(T_vec*1e9, phase, values_from_fit )
            plot2d_2.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')


#    Tabor.set_trigger_source('EVEN')
#    if Pulsing_instrument.get_board_flag():
#        Pulsing_instrument.measurement_close(transfert_info=False)


finally:
    print 'the end'
#    if Pulsing_instrument.get_board_flag():
#        Pulsing_instrument.measurement_close(transfert_info=False)

#    data_measurement.add_data_point(T_vec*1e9, amplitude, phase, real_a, imag_a)
#    data_measurement.close_file()

#    print Pulsing_instrument.measurement_close(transfert_info=True)


#    Tabor.set_trigger_source('EVEN')


#if FIT:
#    data_fit.add_data_point( s.get_fit_params(), s.get_fit_errors())
#    data_fit.close_file()

#plot2d_1.save_png()
#plot2d_2.save_png()
#plot2d_3.save_png()
#plot2d_4.save_png()

#qt.mend()
