from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt
import matplotlib.pyplot as plt
###########################################################
# FIT = False
FIT = True
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
current_source = qt.instruments.get('hp3245')

current_min = -0.5e-3 # in A
current_max = 1.5e-3 # in A +
current_step = 0.05e-3 # in A
# c_ref = 0e-6 # in A
current_vec = np.arange(current_min, current_max+current_step, current_step)

current_source.set_mode('dci')
current_source.set_channel('A')
current_source.set_resolution('high')

SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)
###########################################################

averaging = 500

#StartFrequency [GHz]
f_min_cav =     7.80
#StopFrequency [GHz]
f_max_cav =     7.88
#Step of the Sweep [GHz]
f_step_cav=     0.0002



##############################################################################
power = -0.

freq_vec = np.arange(f_min_cav, f_max_cav + f_step_cav, f_step_cav)
if len(freq_vec) %2 !=0:
        freq_vec = np.arange(f_min_cav, f_max_cav , f_step_cav)

pulse_time = 600.
delta_t = 200.
acq_time = pulse_time + delta_t + 300.

t_rise =None
###########################################################
#
#
#               Experiment
#
#
###########################################################

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(100)
    Pulsing_instrument.write_onetone_pulsessequence( pulse_time*1e-9, t1_start=0.2e-6, m1_start=0.1e-6, delete = 'all',t_rise =t_rise )

qt.mstart()
# qt.msleep(1)
data_measurement = qt.Data(name='Onetone_Spectroscopy_vs_flux')
data_measurement.add_coordinate('Read-out frequency [GHz]', units = 'GHz')
data_measurement.add_coordinate('Current [uA]', units = 'uA')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('S21 normed',            units = 'Volt')
data_measurement.add_value('Phase normed',            units = 'rad')
data_measurement.add_value('Phase unwrap',            units = 'rad')

data_measurement.create_file()

plot3d_1 = qt.Plot3D(data_measurement,
                  name      = 'S21 ',
                  coorddim  = (0,1),
                  valdim    = 2)

# plot3d_2 = qt.Plot3D(data_measurement,
#                     name      = 'Phase ',
#                     coorddim  = (0,1),
#                     valdim    = 3)

# plot3d_3 = qt.Plot3D(data_measurement,
#                   name      = 'S21 normed',
#                   coorddim  = (0,1),
#                   valdim    = 4)
#
# plot3d_4 = qt.Plot3D(data_measurement,
#                     name      = 'Phase normed',
#                     coorddim  = (0,1),
#                     valdim    = 5)

plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 cut',
                  coorddim  = 0,
                  valdim    = 2,
                  maxtraces = 2)

# plot2d_2 = qt.Plot2D(data_measurement,
#                     name      = 'Phase cut',
#                     coorddim  = 0,
#                     valdim    = 3,
#                     maxtraces = 2)

# plot2d_3= qt.Plot2D(data_measurement,
#                     name      = 'Phase unwrap',
#                     coorddim  = 0,
#                     valdim    = 6,
#                     maxtraces = 2)

if FIT:
    fres = np.zeros(len(current_vec))
    f0_err =  np.zeros(len(current_vec))
    data_fit= qt.Data(name='Spectro_fit')
    data_fit.add_coordinate('Current', units = 'uA')
    data_fit.add_value('f_cav ',            units = ' GHz')
    data_fit.add_value('f_cav_error ',            units = 'GHz')
    data_fit.create_file()

    s = fit.S21dB_pic_amplitude()

    # guess parameters##########################################################
    Qi = 1e3
    Qext = 1e3
    f0 = 6.95
    background = 0.
    p0 = [Qi, Qext, f0, background]
    i = 0
i = 0
for c in current_vec:
    current_source.set_current(c)
    qt.msleep(0.1)


    board_flag = None
    try:
        Pulsing_instrument.prep_onetone(freq_vec, averaging, power, acq_time, pulse_time, delta_t)
        qt.msleep(1)
        board_flag = True
        smb_cavity.restartsweep()
        qt.msleep(1)


        Tabor.set_trigger_source('TIM')
        while ats9360.get_completed_acquisition() != 100.:

            result = ats9360.measurement()
            ((real, rea0), (imag,ima0))= result

            real_a = real - rea0
            imag_a = imag - ima0
            real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            amplitude = np.sqrt(real_a**2+imag_a**2)

            complexe = (real_a + 1j*imag_a )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase = np.angle(complexe)

            # (real_a, imag_a) = result
            #
            # real_a = np.mean(np.reshape(real_a, (len(freq_vec),Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            # imag_a = np.mean(np.reshape(imag_a, (len(freq_vec),Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            # amplitude=np.sqrt(real_a**2+imag_a**2)
            #
            #
            # complexe = (real_a + 1j*imag_a )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            # phase = np.angle(complexe)
            phase_unwrap = np.unwrap(phase)
            if i == 0:
                a0 = amplitude
                ph0 = phase
            amp_normed = amplitude - a0
            phase_normed = np.unwrap(phase - ph0)
            qt.msleep(0.1)

            # s21dB = 20*np.log10(amplitude/cos_amplitude_read_out)
            # plot2d_1.replace_inline_data(freq_vec, amplitude)
            # plot2d_2.replace_inline_data(freq_vec, phase)
            # plot2d_3.replace_inline_data(freq_vec, phase_unwrap)

        Tabor.set_trigger_source('EVEN')
        ats9360.measurement_close(transfert_info=False)
        board_flag = False


    finally:
        if board_flag:
            ats9360.measurement_close(transfert_info=False)

        data_measurement.add_data_point(freq_vec,1e6*c*np.ones_like(freq_vec), amplitude, phase, amp_normed,
                phase_normed, phase_unwrap)
        data_measurement.new_block()



        print ats9360.measurement_close(transfert_info=True)
        # smb_cavity.set_freqsweep('OFF')
        # smb_cavity.set_gui_update('ON')
        Tabor.set_trigger_source('EVEN')

    if FIT:
        s.set_data(freq_vec, amplitude)
        if i>0:
            p0 = p
        p = s.fit(p0)
        print 'params:', s.get_fit_params()
        print 'errors:', s.get_fit_errors()
        fres[i] = p[2]
        f0_err[i] = s.get_fit_errors()[2]

        i += 1
        data_fit.add_data_point(c, p[2], s.get_fit_errors()[2])

# hp3245.set_current(c_ref)

data_measurement.close_file()
plot3d_1.save_png()
# plot3d_2.save_png()
qt.mend()

if FIT:

    data_fit.close_file()
    fig, ax =plt.subplots(1,1)
    ax.plot(current_vec*1e6, fres, '+')
    ax.grid()
    ax.set_xlabel('Current [uA]')
    ax.set_ylabel('Resonant_frequency [GHz]')

    plt.show()
