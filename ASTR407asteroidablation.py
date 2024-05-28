#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt


GM = 6.67e-11 * 5.972e24  # Gravitational constant * mass of Earth (m^3/kg/s^2)
Rp = 6371e3  # Radius of Earth (m)
rho0 = 1.225  # Gas density at Earth's Surface (kg/m^3)
T0 = 288.15  # Isothermal Atmosphere termperature (K)
mu = 29e-3  # Molar mass of air (kg/mol)
CD = 1  # Drag coefficient is set to 
rhom = 2500  # Density of the asteroid (kg/m^3)
diameter = 22  # Diameter of the asteroid (m)
radius = diameter / 2  # Radius of the asteroid (m)
CH = 0.1  # Chondrite heat ablation coefficient
Q = 8.26e6  # Heat of ablation (J/kg)
vc = 3000  # Critical velocity (m/s)
v_initial = -20e3  # Initial velocity (m/s)

Rg = 8.314  # Universal gas constant (J/(mol*K))
P0 = Rg * rho0 * T0 / mu  # Surface pressure (Pa)
gp = GM / Rp**2  # Gravitational acceleration at the surface (m/s^2)
H = P0 / (rho0 * gp)  # Scale height (m)

# Functions for atmosphere, drag, gravitational influence, ablation rate
def atmosphere(alt):
    return rho0 * np.exp(-alt / H)

def drag(v, alt, CD, radius):
    rho = atmosphere(alt)
    return -0.5 * CD * rho * np.pi * radius**2 * v * abs(v) / (rhom * (4/3) * np.pi * radius**3)

def grav(alt):
    return -GM / (Rp + alt)**2

#Ablation rate written with respect to the associated heat and the ablation factor which is set to 
# (v^2 - vc^2)/v^2

def ablation_rate(v, alt, radius):
    if abs(v) > vc:
        rho = atmosphere(alt)
        return -CH * rho * v**2 * np.pi * radius**2 / Q
    return 0

# Numerical method that's used to combine the three forces with a time step 0.001.

def Integrate(dt=0.001, radius=radius, vi=v_initial, strength=3e6):
    # First set the velocity to be the incoming velocity and set initial conditions
    v = vi
    alt = 100e3  # Start altitude (100 km)
    mass = rhom * (4/3) * np.pi * radius**3
    alts, vs, radii, energy_deposited = [], [], [], []
    total_energy = 0
    fragmentation_altitude = None  # This variable is used to iterate through altitudes for fragmentation through the integration

    while alt > 0 and mass > 0:
        alts.append(alt)
        vs.append(v)
        radii.append(radius)

        # Change in velocity and corresponding altitude by Newton's method (using derivative)
        rho = atmosphere(alt)
        
        vdot = drag(v, alt, CD, radius) + grav(alt)
        v += vdot * dt
        alt += v * dt

        # Ground impact at altitude = 0 m
        if alt < 0:  
            # Adjust time step to stop exactly at ground level
            dt += alt / v  
            alt = 0

        # Change in mass
        dm = ablation_rate(v, alt, radius) * dt
        mass += dm
        if mass > 0:
            old_radius = radius
            radius = ((3 * mass) / (4 * np.pi * rhom))**(1/3)
            # Kinetic Energy and energy balance
            delta_ke = 0.5 * rhom * (4/3) * np.pi * (old_radius**3 - radius**3) * v**2
            total_energy += delta_ke
            # Convert J to kilotonnes of TNT
            energy_deposited.append(delta_ke / 4.184e12)  

        # Check for fragmentation by dynamical pressure
        dynamic_pressure = 0.5 * rho * v**2
        if dynamic_pressure > strength and fragmentation_altitude is None:
            fragmentation_altitude = alt
        
        if alt <= 0:
            break  # Stop if it reaches or goes below ground level

    max_energy_index = np.argmax(energy_deposited)
    altitude_of_peak_energy_deposition = alts[max_energy_index]
    total_energy_transferred = total_energy / 4.184e12  # Convert J to kilotonnes of TNT

    return alts, vs, radii, energy_deposited, total_energy, altitude_of_peak_energy_deposition, total_energy_transferred, fragmentation_altitude, radii[-1], vs[-1]

# Run Integration and print results
# Store the results of integration in a variable that is updated with the final values
results = Integrate()
alts, vs, radii, energy_deposited, total_energy, altitude_of_peak_energy_deposition, total_energy_transferred, altitude_of_fragmentation, final_radius, final_velocity = results

print(f"The altitude of peak energy deposition: {altitude_of_peak_energy_deposition} m")
print(f"Total energy transferred to the atmosphere: {total_energy_transferred} kilotonnes of TNT")
print(f"The altitude of fragmentation: {altitude_of_fragmentation} m")
print(f"The size of fragments that reach the ground: {final_radius} m")
print(f"The speed of fragments that reach the ground: {final_velocity} m/s")

# Plotting all the results
plt.figure(figsize=(12, 10))

plt.subplot(3, 2, 1)
plt.plot(alts, energy_deposited)
plt.title('Energy Deposition Rate (kilotonnes TNT/s)')
plt.xlabel('Altitude (m)')
plt.ylabel('Energy Rate (kilotonnes TNT/s)')
plt.gca().invert_xaxis()

plt.subplot(3, 2, 2)
plt.plot(alts, np.cumsum(energy_deposited))
plt.title('Cumulative Energy Released (kilotonnes TNT)')
plt.xlabel('Altitude (m)')
plt.ylabel('Cumulative Energy (kilotonnes TNT)')
plt.gca().invert_xaxis()

plt.subplot(3, 2, 3)
plt.plot(alts, radii)
plt.title('Asteroid Size vs Altitude')
plt.xlabel('Altitude (m)')
plt.ylabel('Radius (m)')
plt.gca().invert_xaxis()

plt.subplot(3, 2, 4)
plt.plot(alts, vs)
plt.title('Velocity vs Altitude')
plt.xlabel('Altitude (m)')
plt.ylabel('Velocity (m/s)')
plt.gca().invert_xaxis()

plt.tight_layout()
plt.show()


# 
# ### From the above code and plots, the altitude of peak energy deposition is **3.05 m** with the total energy being transferred is **0.037 kilotonnes of TNT**. Fragmentation occurs as a result of ablation at an altitude of **37.03 km**. I hypothesize that the ablation factor of 0.1 and the heat of about 8260 kJ would lead to a relatively low fragmentation to approximately half sizes. Size of fragments in diameter is about **10.99 m** and the speed at which they reach the ground is **17.4 km/s** which is seen to decrease from the incoming velocity of **20 km/s** as a result of a balance between energy transfer, drag, ablation and gravitational exertion.

# In[ ]:




