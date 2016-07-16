# Import TP profiles from Robinson & Catling, make TP from analytical Sanchez eq ? 
# Save in newTPs folder
# Using cloudBenneke intersection to find transition point

# load things

import numpy as np
import matplotlib.pylab as plt
plt.ion()
from cloudBenneke import cloudBenneke as cb
import loadfile as lf

planets = ['Venus','Earth','Mars','Titan','Jupiter','Saturn','Uranus','Neptune']


def PT(T,P0,T0,g,gamma_a,Rs):
    """                                                                                                           
    Use planet parameters to generate analytic TP profile from Sanchez equation #7                      
    """
    return P0 * (T/T0)**(g/(gamma_a * Rs))


T = np.arange(1,801,5.0)

for planet in planets:
    # Data
    temp,press = lf.getTP(planet)
    mu, Rs, Cp, g, gamma_a,  T0, P0 = lf.PPdat(planet)

    pcross = 0
    tcross = 0
    icross = 0
    for i, p in enumerate(press):
        Ptest = PT(temp[i],P0,T0,g,gamma_a,Rs)
        if np.abs(Ptest - p) < np.abs(pcross - p):
            pcross = Ptest
            tcross = temp[i]
            icross = i
    if planet == 'Mars':
        tcross = temp[-1] + 60.
        pcross = PT(tcross,P0,T0,g,gamma_a,Rs)
        icross = len(temp)

    # Make analytical T,P starting from tcross
    T = np.linspace(tcross,800.0,100)
    P = PT(T,P0,T0,g,gamma_a,Rs)

    # plot figure
    plt.figure()
    plt.semilogy(temp[0:icross],press[0:icross])
    plt.plot(T,P)
    plt.plot(tcross,pcross,'mo')
    plt.ylim(1000,.001)

    # concatenate
    p_all  = np.concatenate((press,P))
    t_all  = np.concatenate((temp,T))

    # write to file
    filename = 'new_PTs/' + planet + '.txt'
    f = open(filename, 'w')
    f.write('# Consolidated Robinson and Catling data files and analytical TP \n')
    f.write('#' + planet + ' \n')
    f.write('# T(K)     Pressure (Bars)  \n')
    for ii in np.arange(0,len(p_all)):
        f.write(np.str(t_all[ii]) + '\t' + np.str(p_all[ii]) + '\n')
    f.close()

