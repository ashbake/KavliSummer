import numpy as np

dirName = '/Users/ashbake/Documents/KavliExoplanetAtm/code/vaporpress/solsys_PTs/'

def PPdat(planet):
    """
    Load planet data
    Returns:
    --------
       data :    
       mu(g/mol) R*(J/g/K) Cp(J/g/K) g(m/s2) gamma_a(K/Km) T0(K) P0(bar)
    """
    f = open('TPparams.txt', 'rU')
    lines = f.readlines()
    f.close()

    parsing = False
    for i in np.arange(len(lines)):
        if lines[i].startswith(planet):
            parsing = True
        else:
            parsing = False
        if parsing:
            data = lines[i].split()
    
    mu, R, Cp, g, gamma_a,  T0, P0 = data[1:len(data)]
    return float(mu),float(R),float(Cp),float(g),float(gamma_a),float(T0),float(P0)


def SPdat(species):
    """
    Read parameters for species
    ('species', 'lnC','L0','Rv','da','db')
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
    return float(lnC), float(L0), float(Rv), float(da), float(db)


def getTP(planet):
    """
    Read TP profiles for each planet:
    Inputs
    ------
         planet:    str
                    Name of the planet to load T, P for
    Returns
    -------
         temp:      array
                    Loaded temperatures in kelvin
         pres:      array
                    Loaded pressures in bars
    """                
    f = open(dirName + planet + '.txt', 'r')
    lines = f.readlines()
    f.close()
    data = lines[1:]
    temp = np.zeros(len(data), dtype=float)
    pres = np.zeros(len(data), dtype=float)
    for j in np.arange(len(data)):
        temp[j] = data[j].split()[1]
        pres[j] = data[j].split()[0]

    return temp,pres



def specInfo(planet):
    """
    Read in table XX from Sanchez with species info
    """
    f = open(dirName + planet + '_species.txt', 'r')
    lines = f.readlines()
    f.close()
    data = lines[1:]
    specs   = np.zeros(len(data), dtype='|S8')
    state   = np.zeros(len(data), dtype='|S8')
    Xc      = np.zeros(len(data), dtype=float)
    Tc      = np.zeros(len(data), dtype=float)
    Pc      = np.zeros(len(data), dtype=float)
    for j in np.arange(len(data)):
        specs[j] = data[j].split()[0]
        state[j] = data[j].split()[1]
        Xc[j]    = data[j].split()[2]
        Tc[j]    = data[j].split()[3]
        Pc[j]    = data[j].split()[4]
    return specs, state, Xc, Tc, Pc


def readSanchez(file, planet, species, state):
    """
    
    Code that reads Sanchez table made from spreadsheet.

    Parameters:
    ------------
    file: string
          Input file that needs to be read

    planet: string
          Planet name as listed in the file
    species: string
          Species name as listed in the file
    state: string
          Species state as listed in the file

    Returns:
    --------
    data: numpy float array
          Array containing all the data available for that planet, species and state

    Modification history:
    ---------------------
    2016-07-09  Jasmina jasmina@nyu.edu Initial implementation
    """
    # open and read whole file
    f = open(file, 'rU')
    lines = f.readlines()
    f.close()

    # start reading planet files for specific mol and state                                    
    parsing = False
    for i in np.arange(len(lines)):
        if lines[i].startswith(planet):
            parsing = True
        elif lines[i].startswith('\t\t\t'):
            parsing = False
        if parsing and lines[i].startswith(species + '\t' + state):
            data = lines[i].split()

    return data
