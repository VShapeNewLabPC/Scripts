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
# smb_atom2       = qt.instruments.get('smb_3')
ats9360        = qt.instruments.get('ats9360')
SSB_cavity      = qt.instruments.get('SSB_cavity')
SSB_atom      = qt.instruments.get('SSB_atom')
SSB_cav2     = qt.instruments.get('SSB_cav2')
RUDAT       = qt.instruments.get('RUDAT_ph')
Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)
#
SSB_atom.set_freq_start(3.5)
SSB_atom.set_freq_stop(7.5)
SSB_atom.set_conversion_loss(6.)
SSB_atom.set_LO_power(15)
SSB_atom.set_band_type(-1)
SSB_atom.set_IF_frequency(0.05)
#
SSB_cav2.set_freq_start(4)
SSB_cav2.set_freq_stop(8)
SSB_cav2.set_conversion_loss(6.)
SSB_cav2.set_LO_power(15)
SSB_cav2.set_band_type(-1)
SSB_cav2.set_IF_frequency(0.0)

averaging = 3e3
power1 = -0. # readout tone
power2 = -0. # spectro tone
power3 = -0. # pi pulse tone

# Frequency RO [GHz]
# f_cav  = 6.912
f_cav  = 7.603

t_protect = 200e-9
t_3 = 6e-6
t_2 = 5e-6
t_2start = 100e-9 + t_3 - t_2
t1_start = 100e-9 + t_3 + t_protect
t_1 = 0.5e-6

delta_t = 200.
acq_time = delta_t + t_1*1e9 + 200

#StartFrequency [GHz]
f_min = 6.1
#StopFrequency [GHz]
f_max = 6.3
#Step of the Sweep [GHz]
f_step = 0.001
# ats9360.set_samplerate(400.)

freq_vec = np.arange(f_min, f_max + f_step, f_step)
if len(freq_vec) %2 !=0:
        freq_vec = np.arange(f_min, f_max , f_step)

att_start = 30.
att_stop = 0.
att_step = -2
Att_vec = np.arange(att_start, att_stop+att_step, att_step)
###########################################################
#
#
#               Experiment
#
#
###########################################################
if Tabor_loading:
    Pulsing_instrument.set_trigger_time(100)
    Pulsing_instrument.write_threetone_pulsessequence( t_3, t_2, t_1, t_2start=t_2start, t_1start = t1_start, delete = 'all')

qt.mstart()

data_measurement = qt.Data(name='3Tones')
data_measurement.add_coordinate('Excitation frequency [GHz]', units = 'GHz')
data_measurement.add_coordinate('Rudat attenuation [dB]', units = 'dB')
data_measurement.add_value('S21 ',            units = 'Volt')
data_measurement.add_value('Phase ',            units = 'rad')
data_measurement.add_value('Re ',            units = 'Volt')
data_measurement.add_value('Im ',            units = 'Volt')

data_measurement.add_value('S21 normed',            units = 'Volt')
data_measurement.add_value('Phase normed ',            units = 'rad')
data_measurement.add_value('Re normed',            units = 'Volt')
data_measurement.add_value('Im normed',            units = 'Volt')
data_measurement.create_file()


plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'S21 ',
                  coorddim  = 0,
                  valdim    = 2,
                  maxtraces = 2)

plot2d_2 = qt.Plot2D(data_measurement,
                    name      = 'Phase ',
                    coorddim  = 0,
                    valdim    = 3,
                    maxtraces = 2)

plot3d_1 = qt.Plot3D(data_measurement,
                  name      = 'S21 vs att',
                  coorddim  = (0,1),
                  valdim    = 2)

plot3d_2 = qt.Plot3D(data_measurement,
                    name      = 'Phase vs att ',
                    coorddim  = (0,1),
                    valdim    = 3)

plot3d_3 = qt.Plot3D(data_measurement,
                  name      = 'S21 normed vs att',
                  coorddim  = (0,1),
                  valdim    = 6)

plot3d_4 = qt.Plot3D(data_measurement,
                    name      = 'Phase normed vs att ',
                    coorddim  = (0,1),
                    valdim    = 7)

#### reference value ###########################################################
board_flag = None
try:
    Pulsing_instrument.prep_threetone(f_cav, f_cav, freq_vec, averaging, power1, power2, power3,
        onesource=1, acq_time=acq_time, pulse_time= t_1*1e9, delta_t=delta_t)
    Tabor.set_ch3_output('OFF')
    qt.msleep(1)
    smb_atom.restartsweep()
    qt.msleep(1)
    board_flag = True
    Tabor.set_trigger_source('TIM')

    while ats9360.get_completed_acquisition() != 100.:
        print ats9360.get_completed_acquisition()

        result = ats9360.measurement()
        ((real, rea0), (imag,ima0))= result

        real_a = real - rea0
        imag_a = imag - ima0
        real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        amplitude = np.sqrt(real_a**2+imag_a**2)

        complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        phase = np.angle(complexe)
        # (real_a, imag_a)= result
        #
        # real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        # imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
        # amplitude = np.sqrt(real_a**2+imag_a**2)
        #
        # complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
        # phase = np.angle(complexe)

        qt.msleep(0.1)

        plot2d_1.replace_inline_data(freq_vec, amplitude)
        plot2d_2.replace_inline_data(freq_vec, phase)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)

    amp0 = amplitude
    phase0 = phase
    real0 = real_a
    imag0 = imag_a

    print ats9360.measurement_close(transfert_info=True)
    Tabor.set_trigger_source('EVEN')
    Tabor.set_ch3_output('ON')


for att in Att_vec:
    RUDAT.set_attenuation(att)
    qt.msleep(1)
    board_flag = None
    try:
        print 'here'
        Pulsing_instrument.prep_threetone(f_cav, f_cav, freq_vec, averaging, power1, power2, power3,
            onesource=1, acq_time=acq_time, pulse_time= t_1*1e9, delta_t=delta_t)
        qt.msleep(1)
        smb_atom.restartsweep()
        qt.msleep(1)
        board_flag = True
        Tabor.set_trigger_source('TIM')

        while ats9360.get_completed_acquisition() != 100.:
            print ats9360.get_completed_acquisition()

            result = ats9360.measurement()
            ((real, rea0), (imag,ima0))= result

            real_a = real - rea0
            imag_a = imag - ima0
            real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            amplitude = np.sqrt(real_a**2+imag_a**2)

            complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase = np.angle(complexe)
            # (real_a, imag_a)= result
            #
            # real_a = np.mean(np.reshape(real_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            # imag_a = np.mean(np.reshape(imag_a, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) ), axis = 1)
            # amplitude = np.sqrt(real_a**2+imag_a**2)
            #
            # complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            # phase = np.angle(complexe)

            qt.msleep(0.1)

            plot2d_1.replace_inline_data(freq_vec, amplitude)
            plot2d_2.replace_inline_data(freq_vec, phase)

        Tabor.set_trigger_source('EVEN')
        ats9360.measurement_close(transfert_info=False)
        board_flag = False


    finally:
        if board_flag:
            ats9360.measurement_close(transfert_info=False)

        data_measurement.add_data_point(freq_vec, att*np.ones_like(freq_vec), amplitude, phase, real_a, imag_a,
                    amplitude - amp0, phase - phase0, real_a - real0, imag_a - imag0)
        data_measurement.new_block()

        print ats9360.measurement_close(transfert_info=True)
        # smb_atom.set_freqsweep('OFF')
        # smb_atom.set_gui_update('ON')
        Tabor.set_trigger_source('EVEN')

data_measurement.close_file()
plot2d_1.save_png()
plot2d_2.save_png()
plot3d_1.save_png()
plot3d_2.save_png()
plot3d_3.save_png()
plot3d_4.save_png()
qt.mend()
