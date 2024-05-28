#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import matplotlib.pyplot as plt


def dm_dr(r, rho_r, dr):
    return 4 * np.pi * r**2 * rho_r * dr #4piρr^2

def dP_dr(rho_r, m, r, dr):
    return -rho_r * m / r**2 * dr

def solve(gamma):
    radius = np.linspace(0, 5, 1000)
    dm = np.zeros_like(radius)
    dP = np.zeros_like(radius)
    m = np.zeros_like(radius)
    P = np.zeros_like(radius)
    rho_r = np.zeros_like(radius)
    P[0] = 1**gamma
    rho_r[0] = 1
    
    for i, r in enumerate(radius):
        
        if i + 1 < len(radius):
            
            r1 = radius[i + 1]

            # calculate change in mass with change in radius
            dm[i + 1] = dm_dr(r1, rho_r[i], r1 - r)
            m[i + 1] = m[i] + dm[i + 1]

            # calculate change in pressure with change in radius
            dP[i + 1] = dP_dr(rho_r[i], m[i + 1], r1, r1 - r)
            P[i + 1] = P[i] + dP[i + 1]

            # calculate new density ratio
            rho_r[i + 1] = P[i + 1]**(1 / gamma)

            # stop if pressure becomes negative
            if P[i] < 0:
                break
    
    return radius[:i], m[:i], P[:i], rho_r[:i], r

gamma_arr = np.array([1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2])
r_max = np.zeros(9)

for i, gamma in enumerate(gamma_arr):
    sol = solve(gamma)
    r_max[i] = sol[4]
    r = sol[0]
    rho_r = sol[3]
    plt.plot(r, rho_r, lw=1, label=f'$\gamma$ = {gamma}')

plt.legend()
plt.xlabel('Radius, r')
plt.ylabel(r'$\rho(r) $- Densities')
plt.title('Plot for ρ(r) - density vs radius for different n, n<=5 and 2>=γ>=1.2')
plt.show()


# In[14]:


plt.scatter(gamma_arr, r_max)
plt.plot(gamma_arr, r_max)
plt.title('Rs for P(r) >0 for different γ')
plt.xlabel('Adiabatic indices, $\gamma$')
plt.ylabel('$R_s$')
plt.show()


# In[10]:


for i, gamma in enumerate(gamma_arr):
    sol = solve(gamma)
    r_max = sol[4]
    r = sol[0]
    rho_r = sol[3]
    plt.plot(r / r_max, rho_r, lw=1, label=f'$\gamma$ = {gamma}')
    
plt.legend()
plt.xlabel('Normalized radius, $r / R_s$')
plt.ylabel(r'$\rho(r) / \rho_c $')
plt.title('Normalized plot ρ(r) vs R with respect to ρ(c) and Rs')
plt.show()


# In[ ]:




