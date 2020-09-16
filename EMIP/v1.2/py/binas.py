## /usr/bin/env python


# ------------------------------------------------
# help
# ------------------------------------------------

"""
Biology, physics, and chemistry constants.
"""


# ------------------------------------------------
# external
# ------------------------------------------------

# ---------------------------------------------------------------
# gonio
# ---------------------------------------------------------------

# DEPRICATED; USE INSTEAD:  math.pi
## defintions for pi :
#pi = 3.1415926535897931

# DEPRICATED; USE INSTEAD:  math.radians()
## factor to convert to radians from degrees:
#deg2rad = pi/180.0     # rad/deg


# ------------------------------------------------
# earth
# ------------------------------------------------

# Radius of earth as used in EMOS library (ECMWF model),
# see for example "jvod2uv.F"
# NOTE: the value 6.375e6 was used in TM #
ae = 6.371e6     # m

# acceleration of gravity:
grav = 9.80665    # m/s2


# ---------------------------------------------------------------
# molecules, mols, etc
# ---------------------------------------------------------------

# Avogadro number
Avog = 6.02205e23      # mlc/mol

#  # Dobson units:
#  real,parameter         ::  Dobs = 2.68668e16    # (mlc/cm2) / DU
  

#
# mole masses of components
# 

# atom masses (table 38)
xm_H     =    1.00790e-3     # kg/mol
xm_N     =   14.00670e-3     # kg/mol
xm_C     =   12.01115e-3     # kg/mol
xm_S     =   32.06400e-3     # kg/mol
xm_O     =   15.99940e-3     # kg/mol
xm_Cl    =   35.45300e-3     # kg/mol
xm_Na    =   22.98980e-3     # kg/mol
xm_Al    =   26.98150e-3     # kg/mol
xm_Si    =   28.08600e-3     # kg/mol
xm_K     =   39.10200e-3     # kg/mol
xm_Ca    =   40.08000e-3     # kg/mol
xm_Mg    =   24.30500e-3     # kg/mol
xm_As    =   47.9e-3         # kg/mol
xm_V     =   50.9e-3         # kg/mol
xm_Cu    =   63.546          # kg/mol
xm_Zn    =   65.4e-3         # kg/mol
xm_Fe    =   55.8e-3         # kg/mol
xm_Ni    =   58.7e-3         # kg/mol
xm_Cd    =  112.40e-3        # kg/mol
xm_Hg    =  200.59e-3        # kg/mol
xm_Pb    =  207.19e-3        # kg/mol

# isotopes:
xm_Rn222 =  222.0e-3         # kg/mol
xm_Pb210 =  210.0e-3         # kg/mol

# molecule masses:
xm_H2O   =  xm_H*2 + xm_O    # kg/mol
xm_O3    =  xm_O * 3         # kg/mol
xm_NO    =  xm_N + xm_O      # kg/mol
xm_NO2   =  xm_N + xm_O * 2  # kg/mol
xm_SO2   =  xm_S + xm_O * 2  # kg/mol
xm_CO2   =  xm_C + xm_O * 2  # kg/mol

# mixtures:
xm_air     =  28.964e-3     # kg/mol : ~80% N2, ~20% O2, ...
xm_seasalt =  74.947e-3     # kg/mol : NaCl, SO4, ..

# sesalt composition:
# (Seinfeld and Pandis, "Atmospheric Chemistry and Physics",
#  table 7.8 "Composition of Sea-Salt", p. 444)
massfrac_cl_in_seasalt  = 0.5504  # (kg Cl )/(kg seasalt)
massfrac_so4_in_seasalt = 0.0768  # (kg SO4)/(kg seasalt)
# other source:
massfrac_na_in_seasalt  = 0.308   # (kg Na )/(kg seasalt)


# ---------------------------------------------------------------
# gas
# ---------------------------------------------------------------
  
# gas constant                        
Rgas = 8.3144     # J/mol/K


# standard pressure
p0 = 1.0e5    # Pa


# ---------------------------------------------------------------
# other
# ---------------------------------------------------------------

# melting point
T0 = 273.16    # K

  
  
  
# ------------------------------------------------
# end
# ------------------------------------------------


