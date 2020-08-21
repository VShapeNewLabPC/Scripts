from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
###########################################################
FIT = False
# FIT = True
Tabor_loading = 0
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
ats9360        = qt.instruments.get('ats9360')
SSB_cavity      = qt.instruments.get('SSB_cavity')
Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(+1)
SSB_cavity.set_IF_frequency(0.0)

averaging = 500
power = -0. #in [dB]
#StartFrequency [GHz]
f_min_cav  = 6.5
#StopFrequency [GHz]
f_max_cav = 7.5


# f_min_cav  = 6.
# f_max_cav = 8.

#Step of the Sweep [GHz]
f_step_cav = 0.01

freq_vec = np.arange(f_min_cav, f_max_cav + f_step_cav, f_step_cav)
if len(freq_vec) %2 !=0:
        freq_vec = np.arange(f_min_cav, f_max_cav , f_step_cav)

pulse_time = 500.
acq_time = 900.
delta_t = 200.
t_rise = None
###########################################################
#
#
#               Experiment
#
#
###########################################################
if Tabor_loading:
    Pulsing_instrument.set_trigger_time(100)
    Pulsing_instrument.write_onetone_pulsessequence( pulse_time*1e-9, t1_start=0.2e-6, m1_start=0.1e-6, delete = 'all', t_rise=t_rise)
    # ats9360.set_acquisition_time(acq_time)

qt.mstart()

data_measurement = qt.Data(name='Spectroscopy_mix')
data_measurement.add_coordinate('Read-out frequency [GHz]', units = 'GHz')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Phase unwrapped',            units = 'rad')

data_measurement.create_file()


plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 onetone',
                  coorddim  = 0,
                  valdim    = 1,
                  maxtraces = 2)

plot2d_2 = qt.Plot2D(data_measurement,
                    name      = 'Phase ',
                    coorddim  = 0,
                    valdim    = 2,
                    maxtraces = 2)


board_flag = None
try:
    Pulsing_instrument.prep_onetone(freq_vec, averaging, power, acq_time, pulse_time, delta_t)
    qt.msleep(1)
    smb_cavity.restartsweep()
    qt.msleep(1)

    board_flag = True

    Tabor.set_trigger_source('TIM')
    while ats9360.get_completed_acquisition() != 100.:
        print ats9360.get_completed_acquisition()

        result = ats9360.measurement()
        ((real, rea0), (imag,ima0))= result

        real_a = real #- rea0
        imag_a = imag #- ima0
        real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        amplitude = np.sqrt(real_a**2+imag_a**2)

        complexe = (real_a + 1j*imag_a )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        phase = np.angle(complexe)
        phase_unwrap = np.unwrap(phase)
        qt.msleep(0.1)

        plot2d_1.replace_inline_data(freq_vec, amplitude)
        plot2d_2.replace_inline_data(freq_vec, phase)
        # plot2d_3.replace_inline_data(freq_vec, phase_unwrap)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)

    data_measurement.add_data_point(freq_vec,amplitude,phase, phase_unwrap)
    data_measurement.close_file()

    print ats9360.measurement_close(transfert_info=True)
    smb_cavity.set_freqsweep('OFF')
    smb_cavity.set_gui_update('ON')
    Tabor.set_trigger_source('EVEN')

if FIT:
    data_fit= qt.Data(name='Spectro_fit')
    data_fit.add_value('parameters ',            units = 'none, none, GHz, Volt')
    data_fit.add_value('errors ',            units = 'none, none, GHz, Volt')
    data_fit.create_file()

    s = fit.Lorentzian()
    s.set_data(freq_vec, amplitude)
    # guess parameters##########################################################
    bg = 0
    area = 1.
    f0 = 7.17
    width = 0.02
    p = [bg, area, f0, width]
    # fitting ##################################################################
    p = s.fit(p)
    values_from_fit = s.func(p)
    print 'params:', s.get_fit_params()
    print 'errors:', s.get_fit_errors()

    data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
    plot2d_1.add(freq_vec, values_from_fit)
    data_fit.close_file()

plot2d_1.save_png()
plot2d_2.save_png()
qt.mend()
