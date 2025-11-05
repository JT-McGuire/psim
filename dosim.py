#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import warnings

warnings.simplefilter('ignore')

subprocess.run(["gcc", "psim.c", "-O3", "-o", "particle_sim", "-I", "/usr/include", "-lm", "-lSDL2", "-lSDL2_ttf"])
time.sleep(0.5)

sim = subprocess.Popen(["./particle_sim"], stdout=subprocess.PIPE)

def to_int(line):
    line = line.split(' ')
    line = line[1:]
    vals = []
    for i in line:
        vals.append(int(i))
    return vals

first_rx = False
last_rx = False

plt.style.use(['dark_background'])
figm, axm = plt.subplots()
axmach = axm.twinx()
figm.suptitle("Relative Channel Flow Measurements")


while True:
    line = sim.stdout.readline()
    if not line:
        break
    
    line = line.decode("utf-8")

    if(line[:5] == "SPEED"):
        speed = to_int(line)
        axm.cla()
        axmach.cla()
        first_rx = True
    elif(line[:5] == "PRESS"):
        pressure = to_int(line)
    elif(line[:5] == "TEMPS"):
        temps = to_int(line)
    elif(line[:5] == "DENSI"):
        density = to_int(line)
    elif(line[:5] == "IGLP:"):
        iglp = to_int(line)
    elif(line[:5] == "AREA:"):
        area = to_int(line)
    elif(line[:5] == "MACH:"):
        mach = to_int(line)
    elif(line[:5] == "ENTRO"):
        entropy = to_int(line)
    elif(line[:5] == "KICK:"):
        print(".", end="", flush=True)
    elif(line[:5] == "SLOWD"):
        print()
        print(line)
        last_rx = True
    
    if(not last_rx and not first_rx):
        plt.pause(0.01)

    if(first_rx and last_rx):
        first_rx = False
        last_rx = False
        axmach.plot(np.array(mach)/1000, color="magenta", lw=3)
        axmach.set_ylabel("Mach Number", color="magenta")
        axmach.tick_params(axis="y", labelcolor="magenta")
        axmach.set_ylim([0, max(mach)/1000*1.1])
        axmach.yaxis.set_label_position("right")
        axmach.tick_params(labelbottom=False)
        axm.plot(np.array(speed), label="X Speed", color="red", lw=2)
        axm.plot(np.array(density), label="Density", color="green", lw=2)
        axm.plot(np.array(pressure), label="Pressure", color="yellow", lw=2)
        axm.plot(np.array(iglp), label="IGL Px", color="orange", lw=2)
        axm.plot(np.array(temps), label="Temperature", color="cyan")
        axm.set_ylim([0, max([max(speed), max(pressure), max(density), max(temps)])*1.1])
        axm.legend()
        axm.tick_params(labelleft=False, labelbottom=False)
        

plt.show()

    
    


    
