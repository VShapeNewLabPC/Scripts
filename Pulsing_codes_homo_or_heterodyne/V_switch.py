import serial
import time

#MAKE A module importable
# Check that it is no double connection to serial!!!
pyserial_open = False

def swit(str):
    try:
        if pyserial_open == False:
            ser = serial.Serial('COM7', 9600, timeout = 500)
        else:
            print 'serial was open'
        pyserial_open = True
    except:
        print 'Error of connection to Arduino'
        return False
    ch = '1'
    # Set state of three switches directional way

    if str == "000":
        ch = 'a'
    elif str == "001":
        ch = 'b'
    elif str == '010':
        ch = 'c'
    elif str == '011':
        ch = 'd'
    elif str == '100':
        ch = 'e'
    elif str == '101':
        ch = 'f'
    elif str == '110':
        ch = 'g'
    elif str == '111':
        ch = 'h'

    #set same stete to all switches
    elif str == 'alloff':
        ch = 'a'
    elif str == 'allon':
        ch = 'h'

    #set one state to one switch
    elif str == '1on':
        ch = '4'
    elif str == '2on':
        ch = '6'
    elif str == '3on':
        ch = '8'
    elif str == '1off':
        ch = '3'
    elif str == '2off':
        ch = '5'
    elif str == '3off':
        ch = '7'

    # if there is no such code:
    else:
        print 'ERROR: Wrong command for switch'
        time.sleep(1)
        ser.close()
        return False
        #return False

    #do the job
    try:
        print 'command is:', ch
        ser.write(ch)
        time.sleep(1)
    finally:
        ser.close()
    return ch
    #return True
