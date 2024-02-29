import numpy as np
import matplotlib.pyplot as plt

# Constants and parameters
mu_dynamic = 0.45  # Dynamic coefficient of friction
k_pad = 1.2  # Thermal conductivity of brake pad (W/mK)
k_disc = 52.0  # Thermal conductivity of brake disc (W/mK)
r1 = 0.05  # Inner disk radius (m)
r3 = 0.1  # Outer radius of the disk (m)
T_ambient = 300.0  # Ambient temperature (K)
initial_velocity = 30.0  # Initial velocity of the vehicle (m/s)
final_velocity = 0.0  # Final velocity of the vehicle (m/s)
braking_time = 10.0  # Braking time (s)
time_steps = 140
radial_steps = 50
pressure_min = 0.8e6  # Minimum pressure (Pa)
pressure_max = 2.5e6  # Maximum pressure (Pa)
youngs_modulus_pad = 1500e6  # Young's modulus of brake pad (Pa)
youngs_modulus_disc = 125e9  # Young's modulus of brake disc (Pa)
poissons_ratio_pad = 0.25  # Poisson's ratio of brake pad
poissons_ratio_disc = 0.28  # Poisson's ratio of brake disc
initial_pad_temp = 300.0  # Initial temperature of brake pad (K)
max_pad_temp = 500.0  # Maximum temperature of brake pad (K)
wear_rate = 2500.0 / 1e3  # Wear rate of brake pad (mm/mile)
normal_load = 10000  # Newtons
sliding_velocity = 2  # meters per second
dt = 0.01  # seconds
dr = 0.001  # meters
k = 60  # Thermal conductivity in W/(m*K), typical for cast iron
cp_disc = 460  # Specific heat capacity in J/(kg*K), approximate for cast iron
rho_disc = 7200  # Density in kg/m^3, typical for cast iron
cooling_coeff = 0.030  # Cooling coefficient, this is a tunable parameter based on experimental data
T_disc = np.full(radial_steps, T_ambient)  # Initial temperature of the disc
T_pad = np.full(radial_steps, initial_pad_temp)  # Initial temperature of the pad
velocity_profile = np.linspace(initial_velocity, final_velocity, time_steps)  # Vehicle velocity profile
contact_area = np.pi * (r3**2 - r1**2)
pressure_min = 0
pressure_max = 1e6

def calculate_hertzian_pressure(r, R, pressure_max):
    return (3 * pressure_max) / (2 * np.pi * (R**2 - r**2))

# Improved friction model based on rate-and-state friction laws
def calculate_friction(v, T_disc, T_pad, sliding_velocity, normal_load):
    a = 1e-2  # State evolution parameter, adjust as needed
    b = 1e-3  # Frictional strength parameter, adjust as needed
    mu_s = 0.8  # Static friction coefficient
    mu_d = 0.6  # Dynamic friction coefficient

    mu = mu_s + (mu_d - mu_s) * np.arctan(v / sliding_velocity)  # Rate-and-state friction law
    mu += a * np.log(1 + (sliding_velocity * np.exp(np.clip((T_disc - T_pad) / b, a_min=-10, a_max=10))) / normal_load)  # State evolution law

    return mu

heat_generation_array = np.zeros(radial_steps)

# Numerical method: Implicit finite difference for heat conduction
def solve_heat_conduction(T, k, cp, rho, heat_gen, ambient_temp, cooling_coeff, dt, dx):
    alpha = k / (cp * rho)  # Thermal diffusivity
    N = len(T)
    A = np.zeros((N, N))
    b = np.zeros(N)
    
    # Build matrix A for the linear system
    for i in range(1, N - 1):
        A[i, i - 1] = A[i, i + 1] = alpha / dx**2
        A[i, i] = -2 * alpha / dx**2 - cooling_coeff
    
    A[0, 0] = A[N - 1, N - 1] = 1  # Boundary conditions
    
    # Update b based on heat generation and cooling
    for i in range(1, N - 1):
        b[i] = -heat_gen[i] - cooling_coeff * (T[i] - ambient_temp)
    
    # Boundary conditions
    b[0] = T[0]
    b[N - 1] = T[N - 1]
    
    # Solve the linear system
    T_new = np.linalg.solve(A, b)
    
    return T_new

# Initialize an array to store the specific temperatures you're interested in
specific_temperatures = []

# Main simulation loop
for t in range(time_steps):
    pressure = pressure_min + (pressure_max - pressure_min) * t / time_steps
    v = velocity_profile[t]
    
    heat_generation_array = np.zeros(radial_steps)

    for r in range(1, radial_steps - 1):
        R = (r3 - r1) / 2

        # Calculate non-uniform pressure distribution using Hertzian contact theory
        p = calculate_hertzian_pressure(r, R, pressure)

        # Calculate sliding friction
        friction = calculate_friction(v, T_disc[r], T_pad[r], sliding_velocity, normal_load)

        # Update heat generation term based on pressure and friction
        heat_generation = (friction * contact_area) * dt
        
        cooling_effect = cooling_coeff * (T_disc[r] - T_ambient)
        net_heat_generation = heat_generation - cooling_effect * contact_area * dt

        # Implicit finite difference update for disc temperature
        T_disc = solve_heat_conduction(T_disc, k_disc, cp_disc, rho_disc, heat_generation_array, T_ambient, cooling_coeff, dt, dr)
        
     specific_temp = T_disc[radial_steps//2]  # Ensure this aligns with your simulation's logic
     specific_temperatures.append(specific_temp)  # Collect this temperature for plotting
        
        # Update pad temperature considering temperature-dependent properties
        #alpha_pad = youngs_modulus_pad / (2 * (1 + poissons_ratio_pad))
        #dT_pad = (v / r) * (T_disc[r] - T_pad[r]) * dt
        #T_pad[r] += dT_pad
        #youngs_modulus_eff = youngs_modulus_pad * (1 - poissons_ratio_pad**2) / (1 + poissons_ratio_pad)
        #k_eff = k_pad * (1 - poissons_ratio_pad)
        #epsilon = 1e-10
        #heat_conduction = k_eff * (T_pad[r + 1] - 2 * T_pad[r] + T_pad[r - 1]) / (dr**2 + epsilon)
        #strain_rate = (v / r) * (1 + poissons_ratio_pad)
        #heat_generation_pad = alpha_pad * youngs_modulus_eff * strain_rate * dT_pad
        #T_pad[r] += heat_conduction + heat_generation_pad
        
# Assuming pressure is calculated in Pascals within your loop
pressures_pa = np.linspace(pressure_min, pressure_max, time_steps)

# Plotting (unchanged)
r_values = np.linspace(r1, r3, radial_steps)
plt.figure(figsize=(10, 6))
plt.plot(pressures_pa, specific_temperatures, '-o', label='Specific Temp vs. Pressure')
plt.xlabel('Pressure')
plt.ylabel('Temperature (K)')
plt.title('Improved Temperature Distribution in Brake Rotor with Sliding Contact and Boundary Conditions')
plt.legend()
plt.show()