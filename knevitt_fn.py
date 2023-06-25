import numpy as np
import pandas as pd
import consts
import multiprocessing
from multiprocessing import Pool
import timeit
import phys_fn as pf
"""
functions for implementing the Knevitt 2014 treatment of 
LMXB outbursts

Knevitt 2014
https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.3087K/abstract

Marchant 2017
https://www.aanda.org/articles/aa/pdf/2017/08/aa30188-16.pdf
"""

############################################
#      Standard Parameter Calculations     #
############################################

def get_disksize(M1,M2,P):
    """
    inputs:
        M1: primary mass, Msol
        M2: secondary mass, Msol
        P:  orbital period, days

    outputs:
        disk size in cm 
    """

    Ph = 24*P; q = M2/M1
    
    a  = 3.53e10 * M1**.333 * (1+q)**.333 * Ph **.666 #cm
    RL = a*0.46*q**-.666 / ( .6*q**-.666 + np.log(1+q**-.333) )
    RD = .7*RL

    return RD

def get_eddington_luminoisty(M1):
    """
    inputs:
        M1: primary mass, Msol
    outputs:
        Eddington luminoisty, erg/s
    """

    return 1.47e39*M1/10.

def get_eddington_mass_transfer(M1):
    """
    inputs:
        M1: primary mass, Msol
    outputs:
        Eddington mass transfer, Msol/yr
    """

    return 2.6e-7*M1/10.

def get_mass_transfer(P):
    """
    inputs:
        P: orbital period, days
    outputs:
        mass donated to BH, Msol/yr
    """

    Ph = P * 24.
    Mdonated        = 1e-10*(Ph/2)**-.666          # Msol/yr
    Mdonated[Ph>=2] = 6e-10*(Ph[Ph>=2]/3)**(5./3.) # Msol/yr

    return Mdonated

def get_diskmass(M1,disksize):
    """
    inputs:
        M1:       primary mass, Msol
        disksize: cm
    outputs:
        maximum disk mass: solar masses
    """

    ac    = 0.01
    mdisk = 2.4e21*ac**-0.85*M1**-0.35*(disksize/1e10)**3.05 #g

    return mdisk / consts.Msol

def get_peak_Mdot(disksize):
    """
    inputs:
        disksize in cm
    outputs:
        peak mass transfer to BH in Msol/yr
    """

    mdp = consts.rho*consts.nu*disksize

    return mdp * consts.YtS/consts.Msol

def get_luminosity(eta,disksize):
    """
    inputs:
        eta
        disksize: cm
    outputs:
        outburst luminosity in erg/s
    """

    L = eta*(consts.c*100)**2*consts.rho*consts.nu*disksize
    
    return L

############################################
#       Outburst Treatment Functions      #
############################################

def get_eta(L,Ledd,Md,Mdedd,f,version):
    """
    inputs:
        L:       luminosity in erg/s
        Ledd:    eddington luminoisty in erg/s 

        Md:      mass transfer rate Msol/yr
        Mdedd:   eddington mass transfer Msol/yr
        version:  type of decline in RIA regime (sharp or gradual)
    outputs:
        eta
    """

    C   = L<=f*Ledd

    if type(Md) == type(1.0):
        Md = Md*np.ones_like(L)
    if type(Mdedd) == type(1.0):
        Mdedd = Mdedd*np.ones_like(L)

    eta    = 0.1*np.ones_like(L)
    if version == 'gradual':
        eta[C] = 0.1*Md[C]/(f*Mdedd[C])
    elif version == 'sharp':
        eta[C] = 0.

    return eta
     
def get_peak_RIA_luminosity(disksize, Ledd, Mdedd, f, version):
    """
    Calculate radiatively inefficient accretion (RIA) luminosity
    
    inputs:
        disksize: cm
        Ledd:     eddington luminoisty in erg/s 
        Mdedd:    eddington mass transfer Msol/yr
        f:        unitless fraction for invoking RIA
        version:  type of decline in RIA regime (sharp or gradual)
    outputs:
        L_RIA: erg/s
    """

    L  = get_luminosity(0.1,disksize)
    Md = get_peak_Mdot(disksize)
    eta = get_eta(L, Ledd, Md, Mdedd,f,version)

    return get_luminosity(eta,disksize)/consts.fcor

def get_outburst_time_vec(disksize,Ledd,Mdedd,f,version,nproc,Lthresh):
    """
    Calculate radiatively inefficient accretion (RIA)
    outburst times. Parallel process.

    vectorized version

    inputs:
        disksize: cm
        Ledd:     eddington luminoisty in erg/s 
        Mdedd:    eddington mass transfer Msol/yr
        f:        unitless fraction for invoking RIA
        nproc:    number of processors
        Lthresh:  threshold luminoisty (erg/s)
    outputs:
        t:        time array, shape ( len(disksize), nsearch )
        l:        luminosity array, shape ( len(disksize), nsearch ) 
        tdet:     detectable time, s
        tout:     outburst time, s
    """
    if type(Ledd) == type(1.0):
        Ledd = Ledd * np.ones_like(disksize)
    if type(Mdedd) == type(1.0):
        Mdedd = Mdedd * np.ones_like(disksize)
    if type(Lthresh) == type(1.0):
        Lthresh = Lthresh * np.ones_like(disksize)

    disksizes = np.array_split(disksize, nproc)
    Ledds     = np.array_split(Ledd,     nproc) 
    Mdedds    = np.array_split(Mdedd,    nproc) 
    Lthreshs  = np.array_split(Lthresh,  nproc) 

    start_time = timeit.default_timer()

    if nproc == 1:
        t,l,tdet,tout,lpeak,job_num = calc_outburst_time_vec([
            1,disksize,Ledd,Mdedd,f,version,Lthresh
        ])

        return t,l,tdet,tout,lpeak,job_num

    else:
        inputs = [ [i,d,l,m,f,version,lt] for i,d,l,m,lt in zip( range(len(disksizes)), disksizes,Ledds,Mdedds,Lthreshs) ]
        pp  = Pool(nproc)
        out = pp.map(calc_outburst_time_vec, inputs)
        pp.close()

        t     = np.concatenate([a[0] for a in out])
        l     = np.concatenate([a[1] for a in out])
        tdet  = np.concatenate([a[2] for a in out])
        tout  = np.concatenate([a[3] for a in out])
        lpeak = np.concatenate([a[4] for a in out])
        nums  = np.array([a[5] for a in out])

    # print("Knevitt Outburst Time, f={}, version={}:".format(f,version),timeit.default_timer() - start_time)

        return t,l,tdet,tout,lpeak,nums

def calc_outburst_time_vec(work):
    """
    Calculate radiatively inefficient accretion (RIA)
    outburst times

    vectorized verion

    inputs:
        work = [
        i: job index
        disksize: cm
        Ledd:     eddington luminoisty in erg/s 
        Mdedd:    eddington mass transfer Msol/yr
        f:        unitless fraction for invoking RIA
        version:  type of decline in RIA regime (sharp or gradual)
        Lthresh:  threshold luminoisty (erg/s)
        ]
    outputs:
        t:        time array, shape ( len(disksize), nsearch )
        l:        luminosity array, shape ( len(disksize), nsearch ) 
        tdet:     detectable time, s
        tout:     outburst time, s
        Lp:       peak luminosity
        job_num
    """

    job_num,disksize,Ledd,Mdedd,f,version,Lthresh = work

    # time for irradiated radius to drop below disk radius
    T  = ((disksize**2)/(3*consts.nu))*np.log(consts.Bm*consts.nu*consts.rho/disksize)

    # correct for long period binaries
    T[T<0] = 0

    # irradiated mass at T
    Mh = ((consts.rho*disksize**3)/3) *np.exp(-3*consts.nu*T/disksize**2)

    # total outburst time
    A1   = (3*consts.nu/consts.Bm)**.5
    A2   = Mh**.5
    tout = T + A2/A1

    # minimum time for search
    tmin = -5*np.ones_like(tout)

    # array for detectable outburst time search
    t = 10**np.linspace(tmin,np.log10(tout),3000)[:-1].T #s

    # generate light curve
    eta   = 0.1*np.ones_like(disksize)
    Lpeak = get_luminosity(eta,disksize)
    tau   = (disksize**2)/(3*consts.nu)

    # save 1d for later
    Lpeak1d   = Lpeak.copy()
    Lthresh1d = Lthresh.copy()

    # reshape
    eta      = np.array([eta]).T
    Lpeak    = np.array([Lpeak]).T
    Mdedd    = np.array([Mdedd]).T
    Ledd     = np.array([Ledd]).T
    Lthresh  = np.array([Lthresh]).T

    Ledd = np.repeat(Ledd, 2999,axis=1)
    # print(Ledd)
    # print(Ledd.shape)

    tau   = np.array([tau]).T

    T     = np.array([T]).T
    A1    = np.array([A1]).T
    A2    = np.array([A2]).T

    L     = Lpeak*np.exp(-t/tau)
    Ld    = 0.1*(consts.c*100.)**2*A1*( A2 - A1*(t-T) )  
    L[t>T] = Ld[t>T]

    if f > 0:

        eta = get_eta(L,Ledd,L,Ledd,f,version)
        L   = eta*L/0.1

    L = L /consts.fcor

    L[L<=0] = 1 # help the log fn

    diff = np.abs(np.log10(L)-np.log10(Lthresh))
    tdet = []
    for trow, Lrow in zip(t, diff):
        tdet.append( trow[np.argmin(Lrow)]  )
    tdet = np.array(tdet)
    tdet[Lpeak1d < Lthresh1d] = 0.0

    return t, L, tdet/86400, tout/86400,Lpeak1d/consts.fcor,job_num

def get_quiesence_time(M1, P,disksize):
    """
    inputs:
        M1:       primary mass, Msol
        P:        orbital period, days
        disksize: cm
    outputs:
        quiesence time, days
    """
    mdisk = get_diskmass(M1, disksize) # Msol
    mdon  = get_mass_transfer(P)       # Msol/yr
    return 365 * mdisk / mdon # days

############################################
#                Probabilities             #
############################################

def calc_prob_outburst_occ(tc,tsurvey):
    """
    Calculate probability of one outburst
    during survey

    inputs:
        tc:      outburst cycle time, days
        tsurvey: time for survey, days
    outputs:
        probability
    """
    
    P = tsurvey/tc
    P[P>1] = 1

    return P

def calc_prob_see_outburst(tdet,tobs):
    """
    Calculate probability of seeing 
    the outburst

    inputs:
        tdet: observable time (days)
        tobs: fiducial time (days)
    outputs:
        probability
    """
    P = np.ones_like(tdet)
    P[tdet<1] = tdet[tdet<1]/tobs

    return P

def get_detection_prob(tdet,tc,tsurvey,tobs):

    """
    Calculate probability of observing 
    the LMXB

    inputs:
        tdet:    observable time
        tc:      outburst cycle time, days
        tsurvey: time for survey, days
        tobs:    fiducial time (days)
    outputs:
        probability
    """
    P1 = calc_prob_outburst_occ(tc,tsurvey)
    P2 = calc_prob_see_outburst(tdet,tobs)

    return P1*P2


####################################
#            Implement It          #
####################################

def execute_k14(bcm_cut,nproc=1,Fthresh = 10*consts.CRAB*1e-3):
    """
    Wrapper function for implementing
    the Knevitt 2014 methods
    """
 
    # Fthresh = 10*consts.CRAB*1e-3
    bcm_cut.loc[:,"Lthresh"] = 4*np.pi*(bcm_cut['distance']*consts.kpc/1000.)**2*Fthresh
    C = (bcm_cut.DC_f09 <1) & (bcm_cut.RRLO_2>1) # locate those in overflow

    # outburst luminosity
    for f in [0.0, 0.05]:
        for version in ['sharp','gradual']:
            if len(bcm_cut[C].Ledd) > 0:
                bcm_cut.loc[C,'Lx_k14_{}_{}'.format(f,version)] = get_peak_RIA_luminosity(
                                bcm_cut[C].Rmax*100.*consts.Rsol, # disksize cm 
                                bcm_cut[C].Ledd, # erg/s
                                bcm_cut[C].Mdedd_k14, # g/s
                                f, 
                                version
                                )
            else:
                bcm_cut.loc[:,'Lx_k14_{}_{}'.format(f,version)] = np.nan
                
    # outburst time      
    for f in [0.0, 0.05]:
        for version in ['sharp','gradual']:
            if len(bcm_cut[C].Ledd) > 0:
                t,l,td,to,lp,job_num = get_outburst_time_vec(
                                    bcm_cut[C].Rmax*100.*consts.Rsol, # cm
                                    bcm_cut[C].Ledd, 
                                    bcm_cut[C].Mdedd_k14,
                                    f,
                                    version,
                                    nproc,
                                    bcm_cut[C].Lthresh
                                    )
                bcm_cut.loc[C,'td_k14_{}_{}'.format(f,version)] = td
                bcm_cut.loc[C,'to_k14_{}_{}'.format(f,version)] = to
                bcm_cut.loc[C,'DC_k14_{}_{}'.format(f,version)] = td / (to + bcm_cut[C]['tq_f09'])
            else:
                bcm_cut['td_k14_{}_{}'.format(f,version)] = np.nan
                bcm_cut['to_k14_{}_{}'.format(f,version)] = np.nan
                bcm_cut['DC_k14_{}_{}'.format(f,version)] = np.nan

    return bcm_cut