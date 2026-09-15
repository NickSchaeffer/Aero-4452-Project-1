import numpy as np
from Project_Functions import coes2RV
from Project_Functions import tle2coes
mu = 398600
inc, RAAN, ecc, omega, theta, n_rev, date = tle2coes("SDO_TLE.txt")

n = n_rev * 2 * np.pi / 86400 
a = (mu / n**2) ** (1/3)
h = np.sqrt(mu*a*(1-ecc**2))
ecc = 0 

R,V = coes2RV(h,ecc,RAAN,inc,omega,theta,mu)

