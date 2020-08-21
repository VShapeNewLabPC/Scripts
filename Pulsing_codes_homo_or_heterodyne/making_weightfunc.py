from lib.math import fit
import ATS9360.DataTreatment as dt
import numpy as np
import qt

###########################################################
#
#
#               Devices
#
#
###########################################################

Tabor           = qt.instruments.get('Tabor')
smb_cavity      = qt.instruments.get('smb_1')
ats9360         = qt.instruments.get('ats9360')
SSB_cavity      = qt.instruments.get('SSB_cavity')
Pulsing_instrument = qt.instruments.get('Pulsing_instrument')
###########################################################
SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)


###########################################################
Tabor_loading = 1

averaging = 200
t_read = 1000e-9
power = -0.
power2 = -2.3
freq_cav = 7.029
t2 = 30e-9 ### Tpi
t1 = 100e-9
t1_start = 0.2e-6
t1_1start = 0.1e-6
t_between = 200e-9

initial_delay = 170

###########################################################
ats9360.set_trigger_delay(initial_delay)


t1_1 = t_read
t1_2 = t_read
t2_start = t1_1start + t1_1 + t_between
t1_2start = t2_start + t2



###########################################################
#
#
#               Experiment
#
#
###########################################################

if Tabor_loading:
    Pulsing_instrument.set_trigger_time(200.)
    Pulsing_instrument.write_IQ_alwayspi_several_RO(t2, t1_1, t1_2, t2_start=t2_start ,
            t1_1start=t1_1start, t1_2start=t1_2start, delta_m1_start=0.1e-6, delete=False)
    # Pulsing_instrument.write_IQ_alwayspi_rising_several_RO(t2, t1_1, t1_2,
    # t_rise=t_rise, t2_start=t2_start , t1_1start=t1_1start,t1_2start=t1_2start, delta_m1_start=0.1e-6, delete=False)
ats9360.set_acquisition_time(300)
# print ats9360.get_nb_sequence()

qt.mstart()
#
data_measurement = qt.Data(name='time_plot')
data_measurement.add_coordinate('time ',         units = 'us')
data_measurement.add_value('V (t) ',         units = 'Volt')
data_measurement.add_value('V std', units = 'Volts')
data_measurement.add_value('V2 (t) ',         units = 'Volt')
data_measurement.add_value('V2 std', units = 'Volts')
data_measurement.add_value('Amp (t) ',         units = 'Volt')


plot2d_1 = qt.Plot2D(data_measurement,
                  name      = 'v(t) ',
                  coorddim  = 0,
                  valdim    = 1 )
# plot2d_1.set_style('linespoints')
# plot2d_2 = qt.Plot2D(data_measurement,
#                   name      = 'v std ',
#                   coorddim  = 0,
#                   valdim    = 2 )

plot2d_3 = qt.Plot2D(data_measurement,
                  name      = 'v2(t) ',
                  coorddim  = 0,
                  valdim    = 3 )
# plot2d_3.set_style('linespoints')

# plot2d_4 = qt.Plot2D(data_measurement,
#                   name      = 'v2 std ',
#                   coorddim  = 0,
#                   valdim    = 4 )
# plot2d_5 = qt.Plot2D(data_measurement,
#                   name      = 'amp ',
#                   coorddim  = 0,
#                   valdim    = 5 )
board_flag = None
try:

    Pulsing_instrument.prep_timing( freq_cav, averaging, power)
    qt.msleep(2)
    Tabor.set_ch4_output('ON')
    amplitude2 = 10**((power2)/10.)
    Tabor.set_ch4_amplitude(2*amplitude2)
    print 'here', Tabor.get_ch4_output(), Tabor.get_ch4_amplitude()
    data_measurement.create_file()
    board_flag = True
    Tabor.set_trigger_source('TIM')

    while ats9360.get_completed_acquisition() != 100.:
        print 'percentage done:' + str(ats9360.get_completed_acquisition())+'%'

        result = ats9360.measurement()

        ((vmean, vstd), (vmean2, vstd2)) = result

        plot2d_1.replace_inline_data_y2(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean, vmean2)
        # plot2d_1.replace_inline_data(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean)
        # plot2d_1.replace_inline_data(np.arange(len(v1))/ats9360.get_samplerate()*1e3, v1)

        # plot2d_2.replace_inline_data(np.arange(len(vstd))/ats9360.get_samplerate()*1e3, vstd)
        plot2d_3.replace_inline_data(np.arange(len(vmean2))/ats9360.get_samplerate()*1e3, vmean2)
        # plot2d_3.replace_inline_data(np.arange(len(v2))/ats9360.get_samplerate()*1e3, v2)

        # plot2d_4.replace_inline_data(np.arange(len(vstd2))/ats9360.get_samplerate()*1e3, vstd2)

        qt.msleep(0.1)
        amp = np.sqrt((vmean-vmean[-1])**2+(vmean2-vmean2[-1])**2)
        # plot2d_5.replace_inline_data(np.arange(len(vstd2))/ats9360.get_samplerate()*1e3, amp)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)
    data_measurement.add_data_point(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean, vstd, vmean2, vstd2, amp)

    print ats9360.measurement_close(transfert_info=True)
    plot2d_1.add(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean2)
Tabor.set_trigger_source('EVEN')
data_measurement.new_block()
vg1 = vmean
vg2 = vmean2
ats9360.set_trigger_delay(initial_delay+(t1_1+t_between+t2)*1e9)
board_flag = None
try:

    Pulsing_instrument.prep_timing( freq_cav, averaging, power)
    qt.msleep(2)
    Tabor.set_ch4_output('ON')
    amplitude2 = 10**((power2)/10.)
    Tabor.set_ch4_amplitude(2*amplitude2)
    print 'here', Tabor.get_ch4_output(), Tabor.get_ch4_amplitude()

    board_flag = True
    Tabor.set_trigger_source('TIM')

    while ats9360.get_completed_acquisition() != 100.:
        print 'percentage done:' + str(ats9360.get_completed_acquisition())+'%'

        result = ats9360.measurement()

        ((vmean, vstd), (vmean2, vstd2)) = result

        plot2d_1.replace_inline_data_y2(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean, vmean2)
        # plot2d_1.replace_inline_data(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean)
        # plot2d_1.replace_inline_data(np.arange(len(v1))/ats9360.get_samplerate()*1e3, v1)

        # plot2d_2.replace_inline_data(np.arange(len(vstd))/ats9360.get_samplerate()*1e3, vstd)
        plot2d_3.replace_inline_data(np.arange(len(vmean2))/ats9360.get_samplerate()*1e3, vmean2)
        # plot2d_3.replace_inline_data(np.arange(len(v2))/ats9360.get_samplerate()*1e3, v2)

        # plot2d_4.replace_inline_data(np.arange(len(vstd2))/ats9360.get_samplerate()*1e3, vstd2)

        qt.msleep(0.1)
        amp = np.sqrt((vmean-vmean[-1])**2+(vmean2-vmean2[-1])**2)
        # plot2d_5.replace_inline_data(np.arange(len(vstd2))/ats9360.get_samplerate()*1e3, amp)

    Tabor.set_trigger_source('EVEN')
    ats9360.measurement_close(transfert_info=False)
    board_flag = False


finally:
    if board_flag:
        ats9360.measurement_close(transfert_info=False)
    data_measurement.add_data_point(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean, vstd, vmean2, vstd2, amp)

    print ats9360.measurement_close(transfert_info=True)
    plot2d_1.add(np.arange(len(vmean))/ats9360.get_samplerate()*1e3, vmean2)
Tabor.set_trigger_source('EVEN')
ve1 = vmean
ve2 = vmean2
data_measurement.close_file()
weightfunc = np.sqrt((ve1-vg1)**2+(ve2-vg2)**2)
ats9360.set_trigger_delay(initial_delay)
qt.mend()
