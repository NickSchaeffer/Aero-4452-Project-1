import numpy as np
from Project_Functions import coes2RV
from Project_Functions import tle2coes
from numpy.linalg import norm
mu = 398600
inc, RAAN, ecc, omega, theta, n_rev, date = tle2coes("SDO_TLE.txt")

n = n_rev * 2 * np.pi / 86400 
a = (mu / n**2) ** (1/3)
h = np.sqrt(mu*a*(1-ecc**2))
ecc = 0 

R,V = coes2RV(h,ecc,RAAN,inc,omega,np.deg2rad(theta),mu)

dr_0 = [0,-100,0]
dv_0 = [0,0,0]

# 1st: hop 100-40km 1 day

# 2nd: football 40-20 km  2 day 

# 3rd: hop 40-1 km 3 days 

# 4th: hop 1km-300m 4 days

# 5th: hop 300m-20m  5 days 

# 6th: oscilation at 20m  6 d

# 7th: vbar approach to 0m 7 day 
