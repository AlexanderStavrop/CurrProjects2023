import matplotlib.pyplot as plt

RES_RT = 100800 

# Open the file
with open('Sensors/Lab_1/Code/thermistor_500.csv', 'r') as f:
    lines = f.readlines()

    # Create lists to store the data
    vo = []
    res = []
    temp = []
    delta_res = []
    
    for line in lines:
        pass
        splitted_line = line.split(';')

        vo.append(float(splitted_line[1]))
        res.append(float(splitted_line[2]))
        temp.append(float(splitted_line[3]))
        delta_res.append(RES_RT - float(splitted_line[2])) 


# Plot the V_out vs Temperature
fig = plt.figure()
plt.plot(vo, temp)
plt.xlabel(r'$V_{out} (V)$')
plt.ylabel('Temperature (C)')
plt.title(r'$V_{out}$ vs Temperature')
plt.savefig("Sensors/Lab_1/Review/arduino_Vo_Temp.eps", format='eps') 
plt.savefig("Sensors/Lab_1/Review/arduino_Vo_Temp.png", format='png') 


# Plot the V_out vs Delta_Resistance
fig = plt.figure()
plt.plot(vo, delta_res)
plt.xlabel(r'$V_{out}$ (V)')
plt.ylabel(r'$\Delta R_T$ (C)')
plt.title(r'$V_{out}$ vs $\Delta R_T$')
plt.savefig("Sensors/Lab_1/Review/arduino_Vo_DeltaTemp.png", format='png')
plt.savefig("Sensors/Lab_1/Review/arduino_Vo_DeltaTemp.eps", format='eps')