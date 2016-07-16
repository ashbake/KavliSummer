import sys, os
import numpy as np
import matplotlib.pyplot as plt
import loadfile as lf
plt.ion()

import cloudBenneke
cb = cloudBenneke.cloudBenneke()



# ::::::::::::
# Define Class 
# ::::::::::::
def VP(T,lnC,Rv,L0, da, db, Xc):
    """
    Claudius Clapeyron Equation curves to match Seager Figure 4.11 
    """
    return np.exp(lnC + (1/Rv)*(-L0/T + da*np.log(T) + db*T))/Xc

def PT(T,P0,T0,g,gamma_a,Rs):
    """
    Use planet parameters to generate analytic TP profile from Sanchez equation #7 
    """
    return P0 * (T/T0)**(g/(gamma_a * Rs))


def alt(planet,temp,press,Rs, Tc, Pc, g ):
    """                                                                                              
    Get the altitude from the pressure and the cloud bottom temp                                     
    """
    # Define p_naught depending on if terrestrial or gas giant                                       
    if planet in ['Earth','Mars','Venus','Titan']:
        p_naught = 1.0   # bar                                                                       
    else:
        p_naught = 10.    # bar                                                                       
        
    # Get Pc by finding closest point in Temp array us
        # right now hard coded by finding intersection or using
        # Sanchez Table III parameters

    # Calculate Zc corresponding to Pc and return whole array                                        
    Zc = (Rs * Tc / g) * np.log(p_naught/Pc)
    Zs =  (Rs * Tc / g) * np.log(p_naught/press)
    return Zc, Zs, p_naught

def getHc(L0,da,Tc,db,Rv,Cp,g):
    """                                                                                                
    Cloud Scale Height Following Seager eq. 4.45                                                       
    """
    L = L0 + da * Tc + db * Tc**2
    return round(Rv * Tc**2 * Cp / g / L , 4)

def getH(Rs,Tc,g):
    """                                                                                                
    Atm Scale Height. Seager 4.46 with T~Tc                                                            
    """
    return round(Rs * Tc / g, 4)



# ::::::::::::::::::::::::::: 
# Plot TP, VP profiles 
# ::::::::::::::::::::::::::: 
# Set up Figure layout

# Planets to Include
planets = ['Venus','Earth','Mars','Titan','Jupiter','Saturn','Uranus','Neptune']

# Planet specific plotting
xlo_lims = [100, 200, 90, 40, 0, 0, 0,0]
xhi_lims = [800, 300, 240, 200, 400, 400, 400, 400]
ylo_lims = [10.0,1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0]
yhi_lims = [0.001,0.001, 1e-6, 0.001, 0.001,0.001,0.001,0.001]

scaleP = False   # Scale by Pressure
scaleH = True  # Scale by height

# Set up figure layout
fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(18,10))

# Define T array for analytical TP profils
T = np.arange(1,801,5.0)

#Step through each planet
ij = 0                   # use step through planets
for i in range(0,2):
    for j in range(0,4):
        # Pick planet, load some parameters
        planet = planets[ij]
        mu, Rs, Cp, g, gamma_a,  T0, P0 = lf.PPdat(planet)

        # Plot TP profiles for planets, convert data mb to bars
        temp, pres = lf.getTP(planet)
        axs[i,j].semilogy(temp,pres,'k-',lw=2)

        # Set up Scale height axis for plotting clouds
        axtemp = axs[i,j].twinx()

        # Step through each species
        specs, states, Xc, Tcc, Pcc = lf.specInfo(planet)
        for k in np.arange(1): #np.arange(len(specs)): 
            # Define species, states
            species = specs[k]
            state   = states[k] 
            
            # Load all variables for this species & state
            lnC, L0, Rv, da, db = lf.SPdat(species)
            
            # Calculate analytical TP profile and VP profile
            VPcurve = VP(T, lnC, Rv, L0, da, db, Xc[k])
            PTcurve = PT(T,P0,T0,g,gamma_a,Rs)

            # Plot PT profile and VP profile
            axs[i,j].plot(T,PTcurve,'r-',lw=1)
            axs[i,j].plot(T,VPcurve,'--',label = species + ' ' + state)
            
            # Get altitude array
            cb.readVP(species)
            cb.cloudbase(temp,pres,Xc[k],10.0)
            Tc = cb.tbase
            Pc = cb.pbase

            # Because of mars. Set to values found by eye if no pbase
            if Tc  == None:
                Tc = Tcc[k]
                Pc = Pcc[k]

            # Get Scale Height info
            Hc = getHc(L0, da, Tc, db, Rv, Cp, g)
            H = getH(Rs,Tc,g)
            Zc, Zs, p_naught= alt(planet,temp,pres,Rs, Tc, Pc,g)

            # Plot point over saturation eq point
            axs[i,j].plot(Tc,Pc,'ko')

            # Plot twin axes, cloud H, Hc shaded
            axtemp.fill_between(T,Zc,(Zc + Hc),facecolor='r',alpha = .3)
            axtemp.fill_between(T,Zc,(Zc + H),facecolor='b',alpha = .05)

        # Scale y axes either based on P or on H
        if scaleP:
            Zlow = (Rs * Tc / g) * np.log(p_naught/ylo_lims[ij])
            Zhi = (Rs * Tc/ g) * np.log(p_naught/yhi_lims[ij])
            axs[i,j].set_ylim(ylo_lims[ij],yhi_lims[ij])
            axtemp.set_ylim(Zlow,Zhi)
        elif scaleH:
            maxHeight = 250.0  #km
            Phi = p_naught/(np.exp(g * maxHeight/(Rs * Tc)))
            Plo = p_naught
            axs[i,j].set_ylim(Plo,Phi)
            axtemp.set_ylim(0.0,maxHeight)      #max height set here
            
        # Set plot labels, axes
        axs[i,j].set_title(planet)
        axs[i,j].set_xlim(xlo_lims[ij],xhi_lims[ij])
        axs[i,j].legend(loc='best')
 
        # Iterate to next planet
        ij += 1
        if ij > len(planets) - 1:
            break

axs[0,0].set_ylabel('P (bar)') 
axs[1,0].set_ylabel('P (bar)') 
axs[1,0].set_xlabel('T (K)')
axs[1,1].set_xlabel('T (K)')
axs[1,2].set_xlabel('T (K)')
axs[1,3].set_xlabel('T (K)')

plt.savefig('plots/solarsys_VPcurves4.png')



