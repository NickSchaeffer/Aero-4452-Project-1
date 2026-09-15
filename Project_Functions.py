import numpy as np
from numpy.linalg import norm 
import math

def coes2RV(h, e, RAAN, i, omega, theta, mu):
    rp = ((h**2) / mu) * (1 / (1 + e * np.cos(theta))) * \
         (np.cos(theta) * np.array([1, 0, 0]) + np.sin(theta) * np.array([0, 1, 0]))

   
    vp = (mu / h) * ((-np.sin(theta)) * np.array([1, 0, 0]) +
                      (e + np.cos(theta)) * np.array([0, 1, 0]))

    C1 = np.array([
        [np.cos(RAAN), -np.sin(RAAN), 0],
        [np.sin(RAAN),  np.cos(RAAN), 0],
        [0,             0,            1]
    ])
    C2 = np.array([
        [1, 0,          0],
        [0, np.cos(i), -np.sin(i)],
        [0, np.sin(i),  np.cos(i)]
    ])
    C3 = np.array([
        [np.cos(omega), -np.sin(omega), 0],
        [np.sin(omega),  np.cos(omega), 0],
        [0,              0,             1]
    ])

    Ctot = C1 @ C2 @ C3 
    R = Ctot @ rp
    V = Ctot @ vp
    return R, V

def UV_prop(R0, V0, mu, dt):
    r0 = norm(R0)
    v0 = norm(V0)
    vr0 = np.dot(R0, V0) / r0

    alpha = 2/r0 - v0**2/mu       

    x = np.sqrt(mu) * abs(alpha) * dt

    tol = 1e-8
    for _ in range(1000):
        z = alpha * x**2

        if z > 0:
            C = (1 - np.cos(np.sqrt(z))) / z
        elif z < 0:
            C = (np.cosh(np.sqrt(-z)) - 1) / (-z)
        else:
            C = 1/2

        
        if z > 0:
            S = (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z))**3
        elif z < 0:
            S = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z))**3
        else:
            S = 1/6

        F = (r0*vr0/np.sqrt(mu)) * x**2 * C + (1 - alpha*r0) * x**3 * S + r0*x - np.sqrt(mu)*dt
        dFdx = (r0*vr0/np.sqrt(mu)) * x * (1 - z*S) + (1 - alpha*r0) * x**2 * C + r0

        dx = F / dFdx
        x = x - dx
        if abs(dx) < tol:
            break

    a = 1/alpha

    z = alpha * x**2
    if z > 0:
        C = (1 - np.cos(np.sqrt(z))) / z
        S = (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z))**3
    elif z < 0:
        C = (np.cosh(np.sqrt(-z)) - 1) / (-z)
        S = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z))**3
    else:
        C = 1/2
        S = 1/6

    f = 1 - (x**2/r0) * C
    g = dt - (1/np.sqrt(mu)) * x**3 * S

    R = f*R0 + g*V0
    r = norm(R)

    fdot = np.sqrt(mu) / (r * r0) * (z*S - 1) * x
    gdot = 1 - (x**2/r) * C

    V = fdot*R0 + gdot*V0

    return R, V


def CW_solve(n,t,dr0,dv0):
    phi_rr = np.array([
        [4-3*np.cos(n*t),0,0],
        [6*(np.sin(n*t)-n*t),1,0],
        [0,0,np.cos(n*t)]
        ])
    
    phi_rv = np.array([
        [(1/n)*np.sin(n*t),(2/n)*(1-np.cos(n*t)),0],
        [(2/n)*(np.cos(n*t)-1), (1/n)*(4*np.sin(n*t)-(3*n*t)),0],
        [0,0,(1/n)*np.sin(n*t)]
                      ])
    phi_vr = np.array([
        [3*n*np.sin(n*t),0,0],
        [6*n*(np.cos(n*t)-1), 0,0],
        [0,0,-n*np.sin(n*t)]
        ]) 
    phi_vv = np.array([
        [np.cos(n*t),2*np.sin(n*t),0],
        [-2*np.sin(n*t),4*np.cos(n*t)-3,0],
        [0,0,np.cos(n*t)]
        ])
    dr = phi_rr@dr0 + phi_rv@dv0
    dv = phi_vr@dr0 +phi_vv@dv0
    return dr,dv

def mean_to_true_anomaly(M_deg, e, tol=1e-10, max_iter=200):
    M = np.radians(M_deg)
    E = M if e < 0.8 else np.pi

    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        f_prime = 1 - e * np.cos(E)
        dE = f / f_prime
        E -= dE
        if abs(dE) < tol:
            break
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                            np.sqrt(1 - e) * np.cos(E / 2))
    theta_deg = np.degrees(theta) % 360
    return theta_deg

def tle2coes(filename):
    with open(filename, "r") as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    line1 = lines[0]
    line2 = lines[1]

    date = line1[18:32].strip() 
    inc = float(line2[8:16])             
    RAAN = float(line2[17:25])           
    ecc = float("0." + line2[26:33].strip())  
    omega = float(line2[34:42])         
    M = float(line2[43:51])             
    n = float(line2[52:63])              
    theta = mean_to_true_anomaly(M, ecc)

    return inc, RAAN, ecc, omega, theta, n, date