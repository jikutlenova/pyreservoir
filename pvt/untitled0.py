# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 11:20:57 2020

@author: jiku
"""

import os
os.chdir("/Desktop/PEG-Python/PVT")
import numpy as np
from pvt_gas import *

pressure = 2010 # psi
temp = 110 # Fahrenheit
sg = 0.7 # specific gravity
x_h2s = 0.07 # mole fraction of H2S in gas
x_co2 = 0.1 # mole fraction of CO2 in gas


# calculate gas z-factor 
pseudo_rho, z_factor = gas_zfactor(temp, pressure, sg, x_h2s, x_co2)
# calculate gas density
rhogas = gas_density(temp, pressure, sg, x_h2s, x_co2)
# calculate gas viscosity 
viscogas = gas_mu(temp, pressure, sg, x_h2s, x_co2)
mp(temp, pressure, sg, x_h2s, x_co2)