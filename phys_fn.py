import pandas as pd
import numpy as np
import timeit
import consts
import knevitt_fn as kv

###############################################
#    Set of functions for adding physical     #
#    quantities and disk instability phsyics  #
###############################################

# Luminosity Calculations: Marchant 2017  https://www.aanda.org/articles/aa/pdf/2017/08/aa30188-16.pdf


##########################
#  Eddington Luminosity  #
#   and Mass Transfer    #
##########################

def get_eta(current_BH_mass, initial_BH_mass):
    """
    Calcualte eta, based on Marchant 2017

    inputs:
        current_BH_mass: Msol
        initial_BH_mass: Msol
    outputs:
        eta: unitless
    """

    where_lower_masses = current_BH_mass < np.sqrt(6)*initial_BH_mass 

    eta_lower_masses = 1 - np.sqrt(1-(current_BH_mass/(3*initial_BH_mass))**2)
    eta = np.where(where_lower_masses, eta_lower_masses, 0.42)

    return eta

def get_Ledd(M, kstar_2):
    """
    Calculating eddington luminoisty, following 
    Marchant 2017

    inputs:
        M:       Msol
        kstar_2: type of the donor star
    outputs:
        Ledd:    erg/s 
    """

    Ledd = np.ones_like(M)

    where_H  = ~np.isin(kstar_2, consts.isHe)
    where_He =  np.isin(kstar_2, consts.isHe)

    Ledd[where_H]  = 1.47e39 * (M[where_H]  / 10.) *  ( (1. + 1.) / 1.7 ) ** -1
    Ledd[where_He] = 1.47e39 * (M[where_He] / 10.) *  ( (1  + 0.) / 1.7 ) ** -1

    return Ledd

def get_Mdedd(current_BH_mass, initial_BH_mass, kstar_2):
    """
    Calculating eddington masss transfer rate, 
    following Marchant 2017

    inputs:
        current_BH_mass: Msol
        initial_BH_mass: Msol
        kstar_2:         type of the donor star
    outputs:
        Mdedd:           Msol/yr
    """

    eta = get_eta(current_BH_mass, initial_BH_mass)
    
    where_H  = ~np.isin(kstar_2, consts.isHe)
    where_He =  np.isin(kstar_2, consts.isHe)

    Medd = np.ones_like(current_BH_mass)

    Medd[where_H]  = 2.6e-7 * (current_BH_mass[where_H]  / 10.0) * ( (1. + 1.) / 1.7) ** -1. * (eta[where_H]  / 0.1) ** -1.
    Medd[where_He] = 2.6e-7 * (current_BH_mass[where_He] / 10.0) * ( (1. + 0.) / 1.7) ** -1. * (eta[where_He] / 0.1) ** -1.

    return Medd

##########################
#    Peak Luminosity     #
##########################

def outburst_Lx_F08(Ledd, P):
    """
    Calculate outburst luminosity following Fragos 2009

    inputs:
        Ledd: eddington luminoisty, erg/s
        P:    orbital periods, days
    outputs:
        outburst luminosity in Chandra band, erg/s
    """

    nbol = 0.8 
    eps = 0.5 #should this be eta?
    
    return nbol * eps * np.minimum(2. * Ledd, 2. * Ledd * P * 24. / 10.)


##########################
#     Disk Properties    #
##########################

def get_mdisk(df):
    """
    Calculate the disk mass for each LMXB

    Fragos 2009 for H and Lasota 2008 for He

    inputs:
        post-processed bcm table
        (needs Rmin and Rmax)
    outputs:
        Mdisk_max: Msol
    """

    def H_mass(M, R):
        Rcm = R * consts.Rsol * 100.
        return 2 * np.pi * 11.4 * (1 / 1e10) ** 0.8 * (Rcm ** 2.8) / 2.8 / consts.Msol
        
        # return 2 * np.pi * 72.4 * M ** -0.19 * Rcm ** 2.92 / (2.92 * (1e11) ** 0.92)#g/cm^-2

    def He_mass(M, R):
        Rcm = R * consts.Rsol * 100.
        return 2 * np.pi * 277  * M ** -0.25 * (1 / (100*consts.Rsol)) ** 0.96 * (Rcm ** 2.96) / 2.96  / consts.Msol
        # return 2 * np.pi * 277 * M ** -0.25 * Rcm ** 2.96 / (2.96 * (1e10) ** 0.96)#g/cm^-2

    Mdisk_max = np.ones_like(df.Rmin)

    use_H  = ~df['kstar_2'].isin(consts.isHe)
    use_He =  df['kstar_2'].isin(consts.isHe)

    for condition, fn in zip([use_H, use_He],[H_mass,He_mass]):
        Rmin = df[condition].Rmin
        Rmax = df[condition].Rmax
        M    = df[condition].mass_1

        Mdisk_max[condition] = fn(M,Rmax) - fn(M,Rmin)

    Mdisk_max = np.array(Mdisk_max) #/ consts.Msol

    return Mdisk_max

###########################
#  Critical Mass Tranfer  #
###########################

def Mcrit_irr_D199(M1, M2, P):
    """
    https://ui.adsabs.harvard.edu/abs/1999MNRAS.303..139D/abstract

    inputs:
        M1, M2 (Msol))
        P      (days)
    outputs:
        Mcrit ( Msol/yr) 
    """
    C=5.e-4
    mdot = (2.e15 ) * M1 ** (0.5) * M2 ** (-0.2) * (P * 24) ** (1.4) * (C/5.e-4)**-0.5
    return 365 * 24 * 60 * 60  * mdot / consts.Msol 


def Mcrit_irr_D199_gen(M1, R):
    """
    https://ui.adsabs.harvard.edu/abs/1999MNRAS.303..139D/abstract

    inputs:
        M1 (Msol)
        R  (Rsol)
    outputs:
        Mcrit ( Msol/yr) 
    """
    R10 = R * 100 * consts.Rsol / 1e10
    mdot = 1.5e15 * M1 **-0.4 * R10 ** 2.1
    return 365 * 24 * 60 * 60  * mdot / consts.Msol 

def Mcrit_M02_He(M1,R):
    """
    https://ui.adsabs.harvard.edu/abs/2002ApJ...564L..81M/abstract

    inputs:
        M1   (Msol))
        R    (Rsol)
    outputs:
        Mcrit ( Msol/yr) 
    """
    R10 = R * consts.Rsol * 100. /1e10
    mdot = 5.9e16 * M1 ** -0.87 * R10 ** 2.62
    return 365 * 24 * 60 * 60  * mdot / consts.Msol 

###########################
#  Fragos 2009 Methods    #
###########################

def outburst_times_F09(df):

    Mdisk_max = df['mdisk_eff']
    mc        = df['mc_f09']
    acc       = df['acc']

    Tq = Mdisk_max / acc #years
    To = Mdisk_max * (acc **2) / ( acc * (mc **2 - acc **2) ) #years
    
    return Tq*365, To*355 #days


###########################
#  Final Probabilities    #
###########################

def get_prob(df):
    # add probabilities
    M = df.mass_1
    disksize = df.Rmax*consts.Rsol*100
    mdisk = kv.get_diskmass(M,disksize)
    df['tq_f09'] = 365.*mdisk/df.acc #days
    
    # knevitt
    TSURVERY = 15 * 365 # days
    TOBS     = 1 #day
    FTHRESH  = 10*consts.CRAB*1e-3
    for f in [0.0, 0.05]:
        for version in ['sharp','gradual']:
            td = df['td_k14_{}_{}'.format(f,version)]
            to = df['to_k14_{}_{}'.format(f,version)]
            tq = df['tq_f09']
            df.loc[:,'prob_k14_{}_{}'.format(f,version)] = kv.get_detection_prob(td,tq+to,TSURVERY,TOBS)
            df.loc[df.acc>df.mc_f09,'prob_k14_{}_{}'.format(f,version)] = 0
    
    # fragos
    C = df.DC_f09 <1
    
    df['Lx_f09'] = outburst_Lx_F08(df.Ledd,df.porb)
    df['f_f09']   = df['Lx_f09'] / (4*np.pi* (df.distance*consts.kpc/1000.)**2 )
    threshold = FTHRESH*np.sqrt(TOBS/df['to_f09'])
    threshold[df['to_f09'] >= TOBS] = FTHRESH

    df.loc[~C,'detect_f09'] = 0.0
    df.loc[(threshold<df['f_f09'])&C,'detect_f09'] = 1.0
    df.loc[(threshold>df['f_f09'])&C,'detect_f09'] = 0.0

    df.loc[~C,'frac_f09'] = 0.0
    df.loc[ C,'frac_f09'] = TSURVERY / df[C]['tq_f09']
    df.loc[df.frac_f09>1,'frac_f09'] = 1

    df['prob_f09'] = df.frac_f09 * df.detect_f09
    return df