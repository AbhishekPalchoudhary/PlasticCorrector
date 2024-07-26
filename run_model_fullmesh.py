### **Plastic corrector algorithm** ###
#### Author: Abhishek Palchoudhary ####
############### 25/07/24 ##############

from plastic_corrector import loading_function, visco_plast
import numpy as np

########################################################################################
# Change load function here :

amplitude=0.47        #peak value of f(t)
n_cycles=2           #number of cycles for loading, change to 1/4 for monotonic loading
num_time_steps=100  #number of steps for time-integration
end_time=1           #upper limit of time history

########################################################################################
#Change elasticity parameters here
E = 75500           #Young's modulus
nu = 0.3             #Poisson ratio
sig_y = 170.         #Yield stress

########################################################################################
#Change non-linear hardening parameters here
b = 19.              #Isotropic parameter
Q = 20.             #Isotropic parameter
C = 127499.           #Kinematic parameter
D = 1334.             #Kinematic parameter

########################################################################################
#Provide von Mises from elastic computation for sig_vm_e = sig_y here
sig_vm_e_fullmesh = np.loadtxt('sigvme_f1.txt')

########################################################################################
# Define batch size for plastic correction based on your RAM limitations
batch_size = 7500




# Get loading function
time, load, peaksvalleys = loading_function(n_cycles, amplitude, num_time_steps, end_time)

# Create an array to hold p_fullmesh
p_fullmesh = np.empty((num_time_steps,len(sig_vm_e_fullmesh)), dtype=sig_vm_e_fullmesh.dtype)

mask = sig_vm_e_fullmesh >= sig_y  # Solution useful for these points, else e_p_n=0
sig_vm_e = sig_vm_e_fullmesh[mask]

n_points = len(sig_vm_e)

print('Number of points for plastic correction',n_points)
print('Number of steps for time integration',num_time_steps)

n_batches = int(np.ceil(n_points / batch_size))

# Initialize arrays for combined results
combined_p = np.zeros((num_time_steps,n_points))

# Process each batch
start_idx = 0
for i in range(n_batches):
    print('Batch no. ', i)
    end_idx = min((i + 1) * batch_size, n_points)

    batch_sig_vm_e = sig_vm_e[start_idx:end_idx]

    batch_n_points = end_idx - start_idx

    # Initialize arrays for this batch
    e_p_0 = np.zeros(batch_n_points)
    s_0 = np.zeros(batch_n_points)
    e_0 = np.zeros(batch_n_points)
    x_0 = np.zeros(batch_n_points)
    p_reinitialize = np.zeros(batch_n_points)

    # Initialize visco_plast object for this batch
    test = visco_plast()

    test.initialise(batch_n_points, e_p_0, s_0, e_0, E, nu, sig_y, b, Q, C, D, x_0, batch_sig_vm_e, p_reinitialize)

    # Integrate using backward Euler for this batch
    test.integrate_backward_euler(batch_n_points, time, load, peaksvalleys)
    
    # Plot results
    # test.plot(time)
    
    # Extract results for this batch
    batch_p, _ , _ , _ = test.extract_qtys(time)

    # Store results for this batch in the combined arrays
    combined_p[:,start_idx:end_idx] = batch_p

    start_idx = end_idx

print('Number of purely elastic points', len(sig_vm_e_fullmesh[~mask]))

p_fullmesh[:,mask] = combined_p
p_fullmesh[:,~mask] = 0

#Save results
np.savetxt(f"p_plasticcorrector.txt", X=p_fullmesh)
