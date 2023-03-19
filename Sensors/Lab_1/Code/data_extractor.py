import serial
import time

DELAY = 500
SERIAL = '/dev/tty/USB0'
BAUD_RATE = 9600

# Create a file to store the data
f = open("thermistor_" + str(DELAY) + ".csv", "w") # C

# Extracting data from the serial port
while True:
    # Open the serial port
    ser = serial.Serial(SERIAL, BAUD_RATE)
    time.sleep(DELAY/1000)
    
    # Read the data from the serial port
    line = ser.readline().decode('utf-8')

    # Write the data to the file
    f.write(line)
    print(line, end='')
