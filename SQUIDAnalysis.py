# -*- coding: utf-8 -*-
"""
@author: cbrph
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
plt.rcParams.update({'font.size': 40})

Area, AreaUncertainty = 17.655, 0.105696  #mm^2
Thickness, ThicknessUncertainty = 25, 1 #nm
Volume= Thickness*1E-7*Area*0.01 #cm^3
#Propogate uncertainty calculation through quadrature
VolumeError = Volume*np.sqrt((ThicknessUncertainty/Thickness)**2
                                   +(AreaUncertainty/Area)**2)  #cm^3

#Number of high field points to take the diamagnetic background from
nfit = 5

def Load_MPMS(filename):
    #Load a Quantum Design MPMS3.dat file

    with open(filename) as f:
        for i, line in enumerate(f):
            if line.strip() == "[Data]":
                    data_start = i
                    break

    df = pd.read_csv(filename, skiprows=data_start+1)

    x = df["Magnetic Field (Oe)"].to_numpy()/1e4
    dx = np.zeros_like(x)  
    y = df["Moment (emu)"].to_numpy()
    dy = df["M. Std. Err. (emu)"].to_numpy()
    
    # Average measurement temperature
    Temperature = df["Temperature (K)"].mean()

    return x, y, dx, dy, Temperature, df

def Find_Coercivity(x, y, dy, N=1000):
    #Find the coercive fields and monte-carlo derived error

    Crossings = np.where(np.diff(np.sign(y)) != 0)[0]

    Hc = []
    HcErrors = []

    for i in Crossings:
        x1, x2 = x[i], x[i+1]
        y1, y2 = y[i], y[i+1]
        dy1, dy2 = dy[i], dy[i+1]
        
        #Linear interpolation
        H_zero = x1 - y1*(x2-x1)/(y2-y1)
        Hc.append(H_zero)
        
        #Monte-carlo
        H_samples = np.empty(N)
        
        for j in range(N):

            y1_rand = np.random.normal(y1, dy1)
            y2_rand = np.random.normal(y2, dy2)

            H_samples[j] = x1 - y1_rand*(x2-x1)/(y2_rand-y1_rand)
        
        HcErrors.append(np.std(H_samples))
        
    #Keep only the two real coercive crossings
    if len(Hc) >= 2:
        HcValue = (np.abs(Hc[0]) + np.abs(Hc[1])) / 2
    else:
        HcValue = np.nan
        
    HcError = np.sqrt(HcErrors[0]*HcErrors[0] + HcErrors[1]*HcErrors[1])/2    

    return HcValue, HcError


# Folder containing the SQUID data
folder = Path("SQUID_data")

Temperatures = []
Magnetisations = []
MagnetisationErrors = []
Coercivities = []
CoercivityErrors = []

plt.figure(figsize=(8,6))
Markers = ['o','s','^','D','v','<','>','p','*','X',
           'P','h','1','2','3','4','8','|','_','+']

files = list(Path("SQUID_data").glob("*.dat"))

#Read and sort by temperature
file_info = []
for file in files:
    _, _, _, _, T, _ = Load_MPMS(file)
    file_info.append((T, file))

file_info.sort(key=lambda x: x[0])

cmap = plt.get_cmap("summer")
colours = cmap(np.linspace(0, 1, len(file_info)))

#Loop over every .txt file
for i, (Temperature, file) in enumerate(file_info):

    print(f"\nProcessing: {file.name}")

    x, y, dx, dy, _, df = Load_MPMS(file)
    print(f"Number of points: {len(x)}")
    
    #Use the highest positive field points
    x_fit = x[-nfit:]
    y_fit = y[-nfit:]

    #Linear fit: y = m*x + c
    Gradient, Intercept = np.polyfit(x_fit, y_fit, 1)
    
    #Background over the whole field range
    Background = Gradient*x

    #Corrected data
    y_corr = y - Background

    Coefficients, CovariantMatrix = np.polyfit(x_fit, y_corr[-nfit:], 1, cov=True)

    GradientAfterCorrection = Coefficients[0]
    SatMoment = Coefficients[1]

    GradientError = np.sqrt(CovariantMatrix[0,0])
    SatMomentError = np.sqrt(CovariantMatrix[1,1])

    Magnetisation = SatMoment/Volume
    MagnetisationError = SatMomentError/Volume
    Coercivity, CoercivityError = Find_Coercivity(x,y_corr,dy)
    
    #Add results to lists
    Temperatures.append(Temperature)
    Magnetisations.append(Magnetisation)
    MagnetisationErrors.append(MagnetisationError)
    Coercivities.append(Coercivity)
    CoercivityErrors.append(CoercivityError)
    
    #    Plot
    # plt.figure(figsize=(6,4))
    # plt.plot(x, y, 'o', ms=10, label='Raw',color ='red')
    # plt.plot(x, Background+Intercept, 'r-', label='Linear background')
    # plt.xlabel('Field')
    # plt.ylabel('Moment')
    # plt.legend()

    plt.plot(
    x*1000,
    y_corr/Volume,
    linestyle='-',
    linewidth=2.5,

    color=colours[i],

    marker=Markers[i % len(Markers)],
    markersize=8,

    markerfacecolor=colours[i], 
    markeredgecolor='black', #colours[i],
    markeredgewidth=2,

    label=f'{Temperature:.1f} K'
)
    #plt.axhline(0, color='k', ls='--')
    
    # plt.errorbar(
    #     x, y_corr,
    #         xerr=dx,
    #         yerr=dy,
    #         fmt='o',
    #         markersize=3,
    #         color ='red',
    #         capsize=2,
    #         elinewidth=0.8,
    #         label='Corrected data'
    #         )

    print("Magnetisation is", Magnetisation, "and gradient is", GradientAfterCorrection)


plt.xlabel('Field (mT)')
plt.ylabel('Magnetisation (emu/cm$^3$)')
plt.legend(title='Temperature', fontsize=30)
plt.show()


Temperatures = np.array(Temperatures)
Magnetisations = np.array(Magnetisations)
Coercivities = np.array(Coercivities)

OrderingIndex = np.argsort(Temperatures)
Temperatures = Temperatures[OrderingIndex]
Magnetisations = Magnetisations[OrderingIndex]
Coercivities = Coercivities[OrderingIndex]

#Systematic uncertainty from the sample volume
Upper = Magnetisations*Volume/(Volume-VolumeError)
Lower = Magnetisations*Volume/(Volume+VolumeError)

plt.figure(figsize=(8,6))

plt.errorbar(
    Temperatures,
    Magnetisations,
    yerr=MagnetisationErrors,
    fmt='H',
    markersize=25,
    markerfacecolor='orange',
    markeredgecolor='#aa0000',
    markeredgewidth=2,
    capsize=4
)

plt.plot(Temperatures, Upper,
         'k--', linewidth=2)

plt.plot(Temperatures, Lower,
         'k--', linewidth=2)

plt.xlabel("Temperature (K)")
plt.ylabel("Magnetisation (emu/cm$^3$)")
plt.title("Temperature dependence of magnetisation")

plt.tight_layout()
plt.show()

plt.figure(figsize=(8,6))

plt.errorbar(
    Temperatures,
    Coercivities*1000,
    yerr=MagnetisationErrors,
    fmt='o',
    markersize=25,
    markerfacecolor='#bcbd22',
    markeredgecolor='blue',
    markeredgewidth=2,
    capsize=4
)

plt.xlabel("Temperature (K)")
plt.ylabel("Coercivity (mT)")
plt.title("Temperature dependence of coercivity")

plt.tight_layout()
plt.show()



