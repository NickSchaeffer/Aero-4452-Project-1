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
dv_1 = -60*n/(-6*np.pi)
# 2nd: football 40-20 km  2 day 
dv_2 = 20*n
# 3rd: hop 40-1 km 3 days 
dv_3 = -39*n/(-6*np.pi)
# 4th: hop 1km-300m 4 days
dv_4 = -.7*n/(-6*np.pi)
# 5th: hop 300m-20m  5 days 
dv_5 = -.280*n/(-6*np.pi)
# 6th: oscilation at 20m  6 d
dv_6  = 0.0005*n
# 7th: vbar approach to 0m 7 day 
T7   = 86400                     
d7   = 0.020                     
dv_7_radial = 2*n*d7            
dv_7_track  = 2*d7/T7            
dv_7 = dv_7_radial + dv_7_track   

dV = 2*(abs(dv_1) + abs(dv_2) + abs(dv_3) + abs(dv_4) + abs(dv_5) + abs(dv_6)) + abs(dv_7)
print(dV*1e3, "m/s")