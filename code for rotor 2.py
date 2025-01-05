import numpy as np

# Parameters
R = 0.15  # Radius of the brake rotor (m)
L = 0.02  # Thickness of the brake rotor (m)
nr = 20  # Number of radial grid points
ntheta = 30  # Number of angular grid points
nz = 10  # Number of axial grid points
alpha = 1e-5  # Thermal diffusivity (m^2/s)
rho = 7800.0  # Density of steel (kg/m^3)
c = 500.0  # Specific heat capacity of steel (J/kg°C)
h = 150.0  # Convective heat transfer coefficient (W/m^2°C)
T_initial = 20.0  # Initial temperature (°C)
T_ambient = 20.0  # Ambient temperature (°C)
P = 1.1e6  # Pressure applied (Pa)z
contact_fraction = 0.25  # Fraction of rotor surface in contact with pad
t_d = 2.0  # Deceleration time (s)
mu = 0.5 #dynamic coefficient of friction
i_v = 100 # initial velocity of the rotor in m/s

# Define custom profiles
def velocity_profile(t):
    # Example: Custom profile function for velocity
    # Modify this function to create different profiles
    if t < t_d:
        return  i_v* (1 - t / t_d)  # Linear decrease
    else:
        return 0.0  # Stop or other conditions

def force_profile(t):
    # Example: Custom profile function for force
    # Modify this function to create different profiles
    if t < t_d:
        return P * (1 - t / t_d)  # Linear decrease
    else:
        return 0.0  # Stop or other conditions

# Initialize temperature array
T = np.ones((nr, ntheta, nz)) * T_initial

# Spatial discretization steps
dr = R / (nr - 1)
dtheta = 2 * np.pi / ntheta
dz = L / (nz - 1)

# Time step and number of steps
dt = 0.01  # Time step (s)
nt = int(t_d / dt)  # Number of time steps

# Perform time integration
for n in range(nt):
    Tn = T.copy()
    t = n * dt
    # Calculate current velocity and force using custom profiles
    v_current = velocity_profile(t)
    F = force_profile(t)
    
    # Avoid negative values
    if F < 0: F = 0
    if v_current < 0: v_current = 0
    
    # Calculate heat generation rate per unit volume
    total_heat_generation = F *mu*v_current * contact_fraction
    volume = np.pi * R**2 * L
    heat_generation_per_volume = total_heat_generation / volume

    # Debug information
    if n % 100 == 0:  # Print every 100 time steps for debugging
        print(f"Time step {n}, Heat Generation: {total_heat_generation:.2f}, Velocity: {v_current:.2f}, Force: {F:.2f}")

    for i in range(1, nr - 1):
        for j in range(1, ntheta - 1):
            for k in range(1, nz - 1):
                # Calculate Laplacian of temperature
                laplacian_T = (
                    (Tn[i+1, j, k] - 2*Tn[i, j, k] + Tn[i-1, j, k]) / dr**2 +
                    (Tn[i, (j+1) % ntheta, k] - 2*Tn[i, j, k] + Tn[i, (j-1) % ntheta, k]) / (R * dtheta)**2 +
                    (Tn[i, j, k+1] - 2*Tn[i, j, k] + Tn[i, j, k-1]) / dz**2
                )

                # Update temperature using the heat conduction equation
                T[i, j, k] = Tn[i, j, k] + alpha * dt * laplacian_T + (heat_generation_per_volume / (rho * c)) * dt

    # Apply convective cooling only on the exposed surfaces:
    # 1. Outer radial surface (r = R)
    for j in range(ntheta):
        for k in range(nz):
            T[-1, j, k] -= (h * (T[-1, j, k] - T_ambient) * dt) / (rho * c * dr)
    
    # 2. Top axial surface (z = 0) and bottom axial surface (z = L)
    for i in range(nr):
        for j in range(ntheta):
            T[i, j, 0] -= (h * (T[i, j, 0] - T_ambient) * dt) / (rho * c * dz)
            T[i, j, -1] -= (h * (T[i, j, -1] - T_ambient) * dt) / (rho * c * dz)

# Calculate the final average temperature over the entire rotor volume
average_temperature = np.mean(T)

print("\nFinal average temperature over the entire rotor volume after braking is:", average_temperature, "°C")

 Generate plot
radial_dist = np.linspace(0, R, nr)
angular_dist = np.linspace(0, 2 * np.pi, ntheta)
R_grid, Theta_grid = np.meshgrid(radial_dist, angular_dist)
temperature_distribution = T[:, :, nz // 2]  # Slice at the mid-plane

# Convert polar coordinates (R, Theta) to Cartesian for plotting
X = R_grid * np.cos(Theta_grid)
Y = R_grid * np.sin(Theta_grid)

plt.figure(figsize=(8, 6))
plt.contourf(X, Y, temperature_distribution.T, cmap="hot")
plt.colorbar(label="Temperature (°C)")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title("Temperature Distribution Across the Brake Rotor")
plt.grid(True)
plt.axis("equal")  # Ensure equal scaling
plt.show()
