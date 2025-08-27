from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

# Create a GEKKO model
m = GEKKO(remote=False)

# Time steps (steady-state model)
m.time = np.linspace(0, 1, 11)

# Constants
R = 0.4615  # kJ/kg.K for water vapor
cp = 2.0    # kJ/kg.K assumed for superheated steam

# Parameters
P_in = m.Param(value=8000)     # kPa, inlet pressure (8 MPa)
T_in = m.Param(value=773.15)   # K, inlet temperature (500 °C)
P_out = m.Param(value=10)      # kPa, outlet pressure
eta = m.Param(value=0.85)      # Turbine isentropic efficiency

# Variables
T_out_s = m.Var(value=450)     # K, outlet temperature (isentropic)
T_out_actual = m.Var(value=500)# K, actual outlet temp
h_in = m.Var(value=cp * T_in.value)  # kJ/kg, enthalpy in
h_out_s = m.Var()              # kJ/kg, isentropic enthalpy out
h_out_actual = m.Var()         # kJ/kg, actual enthalpy out
W_isentropic = m.Var()         # kJ/kg, isentropic work
W_actual = m.Var()             # kJ/kg, actual turbine work

# Equations
m.Equation(h_out_s == cp * T_out_s)
m.Equation(h_out_actual == cp * T_out_actual)
m.Equation(W_isentropic == h_in - h_out_s)
m.Equation(W_actual == eta * W_isentropic)

# Optional constraints or assumptions
m.Equation(T_out_actual >= 320)  # Must stay above saturation temp

# Solve
m.options.IMODE = 4  # Simulation (steady-state)
m.solve(disp=False)

# Results
print(f"Isentropic outlet temperature: {T_out_s.value[0]:.2f} K")
print(f"Actual outlet temperature:     {T_out_actual.value[0]:.2f} K")
print(f"Isentropic work output:        {W_isentropic.value[0]:.2f} kJ/kg")
print(f"Actual turbine work output:    {W_actual.value[0]:.2f} kJ/kg")

# Plot
plt.figure()
plt.plot(m.time, [W_actual.value]*len(m.time), label='Actual Work (kJ/kg)')
plt.plot(m.time, [W_isentropic.value]*len(m.time), '--', label='Isentropic Work (kJ/kg)')
plt.xlabel('Time')
plt.ylabel('Work Output')
plt.legend()
plt.title('Steam Turbine Work Output')
plt.grid(True)
plt.show()
