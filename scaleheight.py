# ::::::::::::::::::::::
# Cloud Scale Height
# ~~~~~~
# Inputs files:
#     
# ::::::::::::::::::::::

import numpy as np

# Load files containing planet and molecule paramters. List of planets will address
# is given in planet_params

# Load Planet Parameters
PPdat = np.genfromtxt('planetTPparams.txt', dtype=None,
                   names=('body','mu', 'Rstar', 'Cp', 'g', 'gamma_a', 'T0', 'P0'))

# Load Species Parameters
SPdat = np.genfromtxt('VPparams.txt', dtype=None,
                   names=('species', 'lnC','L0','Rv','da','db'))


# Load files containing planet TP profile and species info. Save into dictionary
PIdic = {}
for planet in PPdat['body']:
    # Load data for that planet
#    pPT = np.genfromtxt('solsys_PT/' + planet
#                                             + '.txt',names=('pressure','temperature'))
    p_species = np.genfromtxt('solsys_PT/' + planet + 
                                 '_species.txt', dtype=None,
                                 names=('species','state', 'Xc','Tc'))
    for i in np.arange(len(p_species['species'])):
        PIdic[planet] = {}
        PIdic[planet]['species'] = p_species['species']
        PIdic[planet]['state'] = p_species['state']
        PIdic[planet]['Tc'] = p_species['Tc']



class scaleH():
  """
  Load all paramaters for a certain planet and a particular species that will condense in solid or liquid form, and calculate cloud scale height, Hc, and atm scale height, H
  """
  def __init__(self,planet,species,state):
      # Find indices of planet and molecular species
      iplanet = np.where(PPdat['body'] == planet)[0]
      ispecies = np.where(SPdat['species'] == species)[0] 
      self.g  = PPdat['g'][iplanet]          # Planet gravity m/s2
      self.Rv = SPdat['Rv'][ispecies]        # Species Specific R
      iTc = np.where((PIdic[planet]['species'] == species) & (PIdic[planet]['state'] == state))
      self.Tc = PIdic[planet]['Tc'][iTc]     # Planet Cloud Base Temp
      self.L = SPdat['L0'][ispecies]         # Species Latent Heat (J)
      self.Cp = PPdat['Cp'][iplanet]         # Planet specific heat (J)
      self.Rs = PPdat['Rstar'][iplanet]      # Specific Gas Const
  def Hc(self):
      """
      Cloud Scale Height Following Seager eq. 4.45
      """
      return round(self.Rv * self.Tc**2 * self.Cp / self.g / self.L , 4)
  def H(self):
      """  
      Atm Scale Height. Seager 4.46 with T~Tc
      """
      return round(self.Rs * self.Tc / self.g, 4)



# :::::::::::::::::::::::::::::::::
# Save H and Hc to scaleheights.txt
# :::::::::::::::::::::::::::::::::

filename = 'scaleheights.txt'
f = open(filename, 'w')
f.write('# Planets, H and Hc in km, and Tc in Kelvins. \n ')
# Loop through planets and each species
for planet in PPdat['body']:
    f.write(planet + '\n')
    f.write('species'+ '\t' + 'state' + '\t' + 'Hc' + '\t' + 'H' + '\t' + 'Tc' + '\n')
    for i in np.arange(len(PIdic[planet]['species'])):
        species = PIdic[planet]['species'][i]
        state = PIdic[planet]['state'][i]
        f.write(species + '\t' + state + '\t')
        f.write(np.str(scaleH(planet,species,state).Hc()) + '\t')
        f.write(np.str(scaleH(planet,species,state).H()) + '\t')
        f.write(np.str(scaleH(planet,species,state).Tc[0]) + '\n')
    f.write('\n')

f.close()
