import json

Msol     = 1.989e33 #g
SB_sigma = 5.670374419e-8 #SI
c        = 2.99e8 #m/s
Rsol     = 696.34e6 #m
G        = 6.67408e-11 #SI
Lsol     = 3.839e33 #ers/s
YtS      = 365*24*60*60 #s in a year
CRAB     =  2.4e-8
isHe     = [7,8,9,0]
isDonor    = [0,1,2,3,4,5,6,7,8,9] 
isPrimary  = [14] 

final_kstar1 = [14]
final_kstar2 = [0,1,2,3,4,5,6,7,8,9]

# for knevitt
rho  = 1e-8 # g/cm^3
nu   = 5e14 # cm^2/s
kpc  = 3.086e+21 #cm
fcor = 4.
Bm   = 1e5    # cgs

BSEdict_old = {
            'xi': 1.0, 
            'bhflag': 1, 
            'neta': 0.5, 
            'windflag': 3, 
            'wdflag': 1, 
            'alpha1': 1.0, 
            'pts1': 0.001, 
            'pts3': 0.02, 
            'pts2': 0.01, 
            'epsnov': 0.001, 
            'hewind': 0.5, 
            'ck': 1000, 
            'bwind': 0.0, 
            'lambdaf': 0.5, 
            'mxns': 3.0, 
            'beta': 0.125, 
            'tflag': 1, 
            'acc2': 1.5, 
            'remnantflag': 3, 
            'ceflag': 0, 
            'eddfac': 1.0, 
            'ifflag': 0, 
            'bconst': 3000, 
            'sigma': 265.0, 
            'gamma': -1.0, 
            'pisn': 45.0, 
            'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 
            'bhsigmafrac' : 1.0, 
            'polar_kick_angle' : 90.0, 
            'qcrit_array' : [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0, 1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0],
            'cekickflag' : 2, 
            'cehestarflag' : 0, 
            'cemergeflag' : 0, 
            'ecsn' : 2.5, 
            'ecsn_mlow' : 1.4, 
            'aic' : 1, 
            'ussn' : 0, 
            'sigmadiv' :-20.0, 
            'qcflag' : 2, 
            'eddlimflag' : 0, 
            'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 
            'bhspinflag' : 0, 
            'bhspinmag' : 0.0, 
            'rejuv_fac' : 1.0, 
            'rejuvflag' : 0, 
            'htpmb' : 1, 
            'ST_cr' : 1, 
            'ST_tide' : 0, 
            'bdecayfac' : 1, 
            'rembar_massloss' : 0.5, 
            'kickflag' : 0, 
            'zsun' : 0.014}

# f = open("/projects/b1095/siegelj/projects/14_09_014/make/infile.txt","r")
# BSEdict = {}
# for line in f:
#     if line[0] not in [";","[","\n"]:
#         key, val = line.replace("\n","").replace("'","").split(" = ")
#         if key in BSEdict_old.keys():
#             if type(BSEdict_old[key]) == type(1):
#                 BSEdict[key] = int(val)
#             elif type(BSEdict_old[key]) == type(1.0):
#                 BSEdict[key] = float(val)
#             elif type(BSEdict_old[key]) == type([]):
#                 val = val.replace("2.0/21.0",str(2.0/21.0))
#                 BSEdict[key] = json.loads(val)

InitialBinariesColumns  = [
    'kstar_1', 'kstar_2', 'mass_1', 'mass_2', 'porb', 'ecc', 'metallicity',
    'tphysf', 'mass0_1', 'mass0_2', 'rad_1', 'rad_2', 'lum_1', 'lum_2',
    'massc_1', 'massc_2', 'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1',
    'renv_2', 'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2', 'bacc_1',
    'bacc_2', 'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2', 'tms_1', 'tms_2',
    'bhspin_1', 'bhspin_2', 'tphys', 'binfrac','tbirth', 'bin_num','distance',
    'xGx','yGx','zGx','pop']