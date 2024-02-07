import numpy as np
import matplotlib.pyplot as plt

# Given data
AR = 8.90
e = 0.85
K = 1 / (np.pi * AR * e)
CD0 = 0.025
CL_max = 1.7
g = 9.81
W0 = 450300  # N
Thrust0 = 92300  # N
S = 92.5
rho0 = 1.225
TSFC = 0.85 / 3600

# Constants
ft_to_m = 0.3048
lb_to_N = 4.44822
ISA_temp_offset = 15  # ISA+15

# Altitude range
altitude_ft = np.linspace(0, 40000, 100)
altitude_m = altitude_ft * ft_to_m

# Functions
def rate_of_climb(v, altitude):
    rho = rho0 * ((1 - 0.0065 * altitude / 288.15) ** 4.2561)  # Density at given altitude
    CL = (2 * W0) / (rho * v ** 2 * S)  # Lift coefficient
    CL = np.minimum(CL, CL_max)  # Limit lift coefficient to maximum value
    CD = CD0 + K * CL ** 2  # Drag coefficient
    Thrust = Thrust0 * (rho / rho0)  # Thrust at given altitude
    ROC = (Thrust - CD * 0.5 * rho * v ** 2 * S - W0 * g * np.sin(np.arcsin(CL / CL_max))) / W0  # Rate of climb
    ROC = np.nan_to_num(ROC, nan=0, posinf=0, neginf=0)  # Replace invalid values with 0
    return ROC

# Compute rate of climb
ROC_values = []
v_values = np.linspace(0, 400, 100)
for v in v_values:
    ROC = rate_of_climb(v, altitude_m)
    ROC_values.append(np.max(ROC))  # Store the maximum rate of climb at each velocity

# Plot rate of climb
plt.plot(v_values, ROC_values)
plt.xlabel('Horizontal Velocity (m/s)')
plt.ylabel('Rate of Climb (m/s)')
plt.title('Rate of Climb vs Horizontal Velocity')
plt.grid(True)
plt.show()

# Find the best rate of climb and corresponding horizontal velocity
v_best = v_values[np.argmax(ROC_values)]  # Best horizontal velocity for maximum rate of climb
best_ROC = np.max(ROC_values)  # Maximum rate of climb

# Compute time, distance, fuel, and step climb
altitude_step = 40000 * ft_to_m  # Step climb altitude in meters
time_to_climb = abs(altitude_step / best_ROC)  # Take absolute value
distance_traveled = v_best * time_to_climb
fuel_spent = abs(Thrust0 * TSFC * time_to_climb)  # Take absolute value
step_climb = "Yes" if altitude_step > altitude_m[-1] else "No"

# Print results
print("Results:")
print("Best Rate of Climb (m/s):", best_ROC)
print("Best Horizontal Velocity (m/s):", v_best)
print("Time to Climb (s):", time_to_climb)
print("Distance Traveled (m):", distance_traveled)
print("Fuel Spent (N):", fuel_spent)
print("Step Climb Required:", step_climb)
