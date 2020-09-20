def gas_pseudoprops(temp, pressure, sg, x_h2s, x_co2):
  """
  Calculate Gas Pseudo-critical and Pseudo-reduced Pressure and Temperature
  * Pseudo-critical properties
    For range: 0.57 < sg < 1.68
    (Sutton, 1985)
  * Pseudo-reduced properties
    For range: x_h2s (mol%) < 0.738; x_co2 (mol%) < 0.544; 154 < p (psia) < 7026; 40 < temp (°F) < 300 (error 0.97%)
    (Wichert and Aziz, 1972)
  """
  import numpy as np
#
  if sg > 0.57 and sg < 1.68 and x_h2s < 0.738 and x_co2 < 0.544  and pressure > 154 and pressure < 7026 and temp > 40 and temp < 300:
    temp = temp + 459.67 # convert to Rankine

    # calculate pseudocritical properties (Sutton, valid for 0.57<sg<1.68)
    P_pc = 756.8 - (131.07 * sg) - (3.6 * sg**2)
    T_pc = 169.2 + (349.50 * sg) - (74 * sg**2) # in Rankine

    # calculate adjustment to pseudocritical properties for sour gas (Wiechert-Aziz, valid for x_co2<0.544 and x_h2s<0.738)
    e = (120 * (((x_h2s + x_co2)**0.9) - ((x_h2s + x_co2)**1.6))) + (15 * (x_h2s**0.5 - x_h2s**4))
    T_pc = T_pc - e # corrected T_pc
    P_pc = (P_pc * T_pc) / (T_pc - x_h2s * e * (1-x_h2s))

    # calculate pseudoreduced properties
    P_pr = pressure / P_pc
    T_pr = temp / T_pc
  
  else:
    P_pc, T_pc, P_pr, T_pr = np.nan, np.nan, np.nan, np.nan

  return(P_pc, T_pc, P_pr, T_pr)

def gas_zfactor(temp, pressure, sg, x_h2s, x_co2):
  """
  Calculate Gas Compressibility Factor
  For range: 0.2 < P_pr < 30; 1 < T_pr < 3 (error 0.486%)
  (Dranchuk and Aboukassem, 1975)
  """
  # T_pr : calculated pseudoreduced temperature
  # P_pr : calculated pseudoreduced pressure 
  P_pc, T_pc, P_pr, T_pr=gas_pseudoprops(temp, pressure, sg, x_h2s, x_co2)

  from scipy.optimize import fsolve # non-linear solver
  import numpy as np

  if T_pr > 1 and T_pr < 3 and P_pr > 0.2 and P_pr < 30:
    a1 = 0.3265; a2 = -1.0700; a3 = -0.5339; a4 = 0.01569; a5 = -0.05165; a6 = 0.5475
    a7 = -0.7361; a8 = 0.1844; a9 = 0.1056; a10 = 0.6134; a11 = 0.7210

    def f(y):
      rho_pr, z = y
      c1 = a1 + (a2/T_pr) + (a3/(T_pr**3))+ (a4/(T_pr**4))+ (a5/(T_pr**5))
      c2 = a6 + (a7/T_pr) + (a8/(T_pr**2))
      c3 = a9*((a7/T_pr) + (a8/(T_pr**2)))
      c4 = (a10)*(1+(a11*(rho_pr**2)))*((rho_pr**2)/(T_pr**3))*(np.exp(-a11*(rho_pr**2)))

      f1 = z + (c3*(rho_pr**5)) - (c2*(rho_pr**2)) - (c1*(rho_pr**1)) - c4 - 1
      f2 = rho_pr - ((0.27 * P_pr) / (z * T_pr))
      return[f1, f2]

    pseudo_rho, z_factor = fsolve(f, [1, 1]) # initial guess
  
  else:
    pseudo_rho, z_factor = np.nan, np.nan

  return(pseudo_rho, z_factor) # result is density, z-factor

def gas_density(temp, pressure, sg, x_h2s, x_co2):
  """
  Calculate Gas Density
  For range: this is not a correlation, so valid for infinite intervals
  """ 
  red_den,z=gas_zfactor(temp, pressure, sg, x_h2s, x_co2)
  temp = temp + 459.67
  R = 10.732 # gas constant in (ft3*psi)/(lb-mol*R) 
  rhogas = (28.97 * sg * pressure) / (z * R * temp)
  return rhogas 
 
def gas_mu(temp, pressure, sg, x_h2s, x_co2):
  """
  Calculate Gas Viscosity 
  For gas with CO2 and N2 composition
  For range: 100 < temp (°F) < 340; 0.9 < x_CO2 (mol%) < 3.2; x_N2 (mol%) < 4.8 (std 2.7-9.0%)
  (Lee et al, 1996)
  """
  rhogas=gas_density(temp, pressure, sg, x_h2s, x_co2)
  import numpy as np

  if temp > 100 and temp < 340:
    temp = temp + 459.67
    Mg = 28.97 * sg
    rhogas_lee = rhogas * 0.0160185 # lbm/ft3 converted to gas density unit of Lee et al (g/cm3)
    K = ((0.00094 + 2E-06)*(temp**1.5)) / (209 + 19*Mg + temp)
    x = 3.5 + (986 / temp) + (0.01 * Mg)
    y = 2.4 - 0.2*x  
    viscogas = K * np.exp(x * (rhogas_lee**y))
  
  else:
    viscogas = np.nan
  return viscogas

#calculate Al-Hussainy pseudopressure
def mp(temp, p, sg, x_h2s, x_co2):
    pp = 0
    pold = 0
    Xold = 0
    pstep = p / 20
    for N in range(1,20):
        pnew = pold + pstep
        if pnew<154:
            pold=pnew
            pnew=pold+pstep
            N=N+1
        red_den,z=gas_zfactor(temp, pnew, sg, x_h2s, x_co2)
        visg=gas_mu(temp, pnew, sg, x_h2s, x_co2)
        Xnew = 2 * pnew / z / visg
        pp = pp + (Xold + Xnew) / 2 * pstep
        pold = pnew
        Xold = Xnew
    return pp

def gas_fvf(z, temp, pressure):
  """
  Calculate Gas FVF
  For range: this is not a correlation, so valid for infinite intervals
  """
  temp = temp + 459.67
  Bg = 0.0282793 * z * temp / pressure 
  return(Bg)

def gas_fvf2(unit='unit1', z=0.8, temp=186, pressure=2000):
  """
  Gas FVF calculated in other units
  unit: choice of units (unit1: RB/scf, unit2: res m3/std m3)
  for unit1, inputs temp in Rankine (Fahrenheit + 460), pressure in psia or psig
  for unit2, inputs temp in Kelvin, pressure in psia or psig
  """
  if unit == 'unit1':
    return(0.00503676 * z * temp / pressure) 
  if unit == 'unit2':
    return(0.350958 * z * temp / pressure)



def gas_compressibility(T_pr, P_pr, rho_pr, z, P_pc):
  """
  Calculate Gas Isothermal Compressibility
  For range: unspecified
  (Trube, 1957; Mattar, 1975)
  """
  import numpy as np

  a1 = 0.3265; a2 = -1.0700; a3 = -0.5339; a4 = 0.01569; a5 = -0.05165; a6 = 0.5475
  a7 = -0.7361; a8 = 0.1844; a9 = 0.1056; a10 = 0.6134; a11 = 0.7210

  do = ((a1 + (a2/T_pr) + (a3/T_pr**3) +(a4/T_pr**4) + (a5/T_pr**5)) * rho_pr) + \
      (2 * ((a6 + (a7/T_pr) + (a8/T_pr**2))) * rho_pr**2) - \
      (5 * a9 * (((a7/T_pr) + (a8/T_pr**2))) * rho_pr**4) + (1 + (a11 * rho_pr**2) - (a11 * rho_pr**2)**2) \
      * ((2 * a10 * rho_pr / T_pr**3)*np.exp(-a11 * rho_pr**2))

  c_pr_analytical = (1 / P_pr) - ((0.27 / (z**2 * T_pr)) * (do / (1 + ((rho_pr / z) * do))))
  cgas_analytical = c_pr_analytical / P_pc
  return(cgas_analytical)

          
