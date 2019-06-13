import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

hbarc=0.1973

############################################################
####################### EOS data ###########################
############################################################

# Load EOS from MUSIC (generated with https://github.com/j-f-paquet/eos_maker - SMASH branch )
eos_location="./hrg_hotqcd_eos_binary.dat"

raw=np.fromfile(eos_location, dtype=(float,4))

e=raw[:,0] # GeV/fm^3
p=raw[:,1] # GeV/fm^3
s=raw[:,2] # fm^-3
T=raw[:,3] # GeV

cs2=InterpolatedUnivariateSpline(e,p,k=1).derivative(n=1)(e)

# cs2_fct takes the temperature in fm^-1
cs2_qcd_fct=InterpolatedUnivariateSpline(T/hbarc,cs2)
