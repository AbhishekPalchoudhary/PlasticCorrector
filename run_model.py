### **Plastic corrector algorithm** ###
#### Author: Abhishek Palchoudhary ####
############### 25/07/24 ##############

from plastic_corrector import loading_function, visco_plast
import numpy as np

########################################################################################
# Change load function here :

amplitude=0.8        #peak value of f(t)
n_cycles=2           #number of cycles for loading, change to 1/4 for monotonic loading
num_time_steps=1000  #number of steps for time-integration
end_time=1           #upper limit of time history

########################################################################################
#Change elasticity parameters here
E = 200000           #Young's modulus
nu = 0.3             #Poisson ratio
sig_y = 100.         #Yield stress

########################################################################################
#Change non-linear hardening parameters here
b = 10.              #Isotropic parameter
Q = 100.             #Isotropic parameter
C = 40000.           #Kinematic parameter
D = 400.             #Kinematic parameter

########################################################################################
#Provide von Mises from elastic computation for sig_vm_e = sig_y here
sig_vm_e = np.array([489.249, 698.774, 194.126])

# Get loading function
time, load, peaksvalleys = loading_function(n_cycles, amplitude, num_time_steps, end_time)
# Initialize visco_plast object for this batch
n_points = len(sig_vm_e)
print('Number of points for plastic correction',n_points)
print('Number of steps for time integration',num_time_steps)
e_p_0 = np.zeros(n_points)
s_0 = np.zeros(n_points)
e_0 = np.zeros(n_points)
x_0 = np.zeros(n_points)
p_reinitialize = np.zeros(n_points)
test = visco_plast()
test.initialise(n_points, e_p_0, s_0, e_0, E, nu, sig_y, b, Q, C, D, x_0, sig_vm_e, p_reinitialize)
# Integrate using backward Euler for this batch
test.integrate_backward_euler(n_points, time, load, peaksvalleys)
# Plot results
test.plot(time)
# Extract results for this batch
p,ep,s,e=test.extract_qtys(time)
# Save results
np.savetxt(fname=f"p_plasticcorrector.csv", delimiter=";", X=p)
#########################################################################################
