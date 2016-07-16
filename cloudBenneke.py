# Code for class that will go into haze.py in pyrat-bay
# Questions:
# where in PB will get T, P
# where will get X_c
# how will do parameter loading for each species ... molecule.dat

import numpy as np
import loadfile as lf
import matplotlib.pylab as plt

#T,P = lf.getTP(planet)
#cb = cloudBenneke()
#cb.readVP(species)
#cb.cloudbase(T,P,Xc,10.0)

#plt.semilogy(T,P)
#plt.plot(T,cb.VP)
#plt.plot(cb.tbase,cb.pbase,'mo')
#plt.ylim(1000,0.01)

def extinction(pyrat):
    """
    Calculate the cloud extinction coefficient.
    """
    # Load parameters for that molecular species (in future pyrat.cond.names for condensate species) 
    species = 'H2O'
    pyrat.haze.cloudBenneke.readVP(species)   # do this outside once in beginning of code and call it??
    # Fill in Cloud Base Pressure and Temperature
    pyrat.haze.cloudBenneke.cloudbase(pyrat.atm.temp,pyrat.atm.press)
    # Calculate the extinction coefficient (in cm2 molecule-1)
    pyrat.haze.cloudBenneke.extinction(pyrat.spec.wn, pyrat.atm.press)
    


class cloudBenneke():
    """
    Benneke Model for parameterized cloud
    """
    def __init__(self):
        self.name  = "cloud_benneke"  # Model name
        self.pars  = [ 1.0,           # p_top  : Pressure at top of cloud
                       4.0,           # q_star : Condensate fraction one scale height below top
                       2.0,           # Hc     : Cloud shape parameter
                       2.0 ]          # reff   : Effective Particle radius
        self.npars = len(self.pars)   # Number of model fitting parameters
        self.ec    = None             # Model extinction coefficient
        self.pbase = None             # Cloud base pressure (found in intersection fxn)
        self.tbase = None
        self.lnC   = None
        self.L0    = None
        self.Rv    = None
        self.da    = None
        self.db    = None

    def extinction(self, wn, pressure):
        """
        Parameters
        ----------
        wn:  1D float ndarray
        Wavenumber array in cm-1.
        """
        pass

    def readVP(self,species):
        """
        Read vapor pressure parameters file and fill in parameters
        that are needed for vapor pressure curve
        Defines:
        'species', 'lnC','L0','Rv','da','db'
        """          
        f = open('VPparams.txt', 'rU')
        lines = f.readlines()
        f.close()
      
        parsing = False
        for i in np.arange(len(lines)):
            if lines[i].startswith(species):
                parsing = True
            else:
                parsing = False
            if parsing:
                data = lines[i].split()
                
        lnC, L0, Rv, da, db =  data[1:len(data)]
        self.lnC, self.L0, self.Rv, self.da, self.db = \
            float(lnC), float(L0), float(Rv), float(da), float(db)

    def cloudbase(self,T,P,Xc,refpressure):
        """
        Determine the intersection between PT profile and Vapor Pressure Curves
        Returns:
           - pbase:   
               the cloud base pressure for use in extinction
               None if there is no cloud/intersection
           - VP:
               the vapor pressure curve for that species
        Method Notes:
        determine VP: Claudius Clapeyron Equation from Seager/Sanchez
        ~ VP and PT should be on same temperature grid. Find first point of intersection
        starting from surface (high pressure end). Then interpolate arrays around that 
        point and find intersection
        ~ Need T and P, calculate VP curve using T and parameters imported from readVP
        """
        self.VP = np.exp(self.lnC + (1/self.Rv)*(-self.L0/T + self.da*np.log(T) + self.db*T))/Xc

        # Find places where Vapor pressure becomes higher than Pressure
        crossover = np.pad(np.diff(np.array(self.VP > P).astype(int)), (1,0), \
                               'constant', constant_values = (0,))

        # Put places where P < VP to VP being lower. Find 1's
        ipbases = np.argwhere(crossover == 1).reshape(-1)

        # Find intersection at highest P, break if nonexistent
        if len(ipbases) == 0:
            # Check if base should be at bottom of atm: if VP > P
            if np.array(self.VP < P)[np.argwhere(P == max(P)).reshape(-1)] == True:
                self.pbase = P[-1]
                self.tbase = T[-1] # how to get temp?
                print 'Cloud base is the bottom of the atm'
            else:
                self.pbase = None
                self.tbase = None
            # Return Print statement if this occurs
                print 'No cloud base found'
        else:
            ipbase = ipbases[np.where(P[ipbases] == np.max(P[ipbases]))]

        # Pick out subsection
            subV  = np.log10(self.VP[ipbase - 1 : ipbase + 1])
            subP  = np.log10(P[ipbase - 1: ipbase + 1])
            subT  = T[ipbase - 1: ipbase + 1]
            print subV
            print subP
            print subT
            
        # Find intersection using algebra of intersection of two lines, same x0
            m_V = np.diff(subV) / np.diff(subT)
            m_P = np.diff(subP) / np.diff(subT)
            self.pbase = (10**((subV[0] - subP[0] * m_V / m_P) / (1 - m_V / m_P)))[0]
            self.tbase = ((np.log10(self.pbase) - subP[0]) / m_P + subT[0])[0]
        

      


      



