import serial
import time

DELAY = 500
BAUD_RATE = 9600

# Create a file to store the data
f = open("Photoresistor.csv", "w")

# Extracting data from the serial port
while True:
    # Open the serial port
    ser = serial.Serial('/dev/tty/USB0', 9600)
    time.sleep(DELAY/1000)
    
    # Read the data from the serial port
    line = ser.readline().decode('utf-8')

    # Write the data to the file
    f.write(line)
    print(line, end='')
