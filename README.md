## A plastic correction algorithm for full-field elasto-plastic finite element simulations

This algorithm pertains to the numerical implementation of a new local plastic correction algorithm that is aimed to
accelerate elasto-plastic finite element (FE) simulations for structural problems
exhibiting localised plasticity (around e.g. notches, geometrical defects) due to proportional global loading sequences. 
The proposed method belongs to the category of generalised multi-axial Neuber-type
methods, which process the results of an elastic prediction point-wise in order to
calculate an approximation of the elasto-plastic solution. The publication associated to this algorithm can be found here [[1]](#1).

Usage (correction of few elastic points)

In **run_model.py**:

1. Add the von Mises stress values from an elastic computation $`\bar{\sigma}_{\textrm{VM}}^{\#}`$ (at $`f(t)`$ = 1) in the array _sig_vm_e_. This should correspond to $`\bar{\sigma}_{\textrm{VM}}^{\#} = \sigma_y`$, where $`\sigma_y`$ is the yield stress of the material.
2. Set the form of the proportional load function $`f(t)`$ (_amplitude_, _n_cycles_, _num_time_steps_ and _end_time_). For monotonic loading, n_cycles can be set to 1/4.
3. Set the hardening parameters. The elasticity parameters should match that of the elastic computation used to get $`\bar{\sigma}_{\textrm{VM}}^{\#}`$.
4. The results for the quantity of interest (QoI) will be written for the entire time history. The cumulative plastic strain (p) is written by defaut, the user is free to change the QoI as per need.

Usage (correction of a full elastic FEA computation)

In **run_model_fullmesh.py**:


1. The von Mises stress values from an elastic computation $`\bar{\sigma}_{\textrm{VM}}^{\#}`$ (at $`f(t)`$ = 1) is set in the array _sig_vm_e_ via a text file.
   Steps 2-4 remain the same.
5. **set_values_to_mesh.py** may optionally be used to set the plastic corrected values to the mesh. Fenics [[2]](#2) is required for this step. 


## Example of usage with a mesh of a specimen containing pores

We demonstrate the usage of the plastic corrector with the following example boundary value problem (BVP).

A specimen with a sub-volume of pores explicitly meshed is considered. These pores arise due to a casting manufacturing process, and information on their geometry was obtained from computed tomography.

The material parameters (identified by other authors [[3]](#3)) correspond to an aluminium (AlSi7Mg0.3) alloy:

Elasticity parameters:

$`E =`$ 75500 MPa, 
$`\sigma_y =`$ 170 MPa. 

Isotropic and kinematic hardening parameters:

$`b =`$ 19, 
$`Q =`$ 20 MPa, 
$`C =`$ 127499 MPa, 
$`D =`$ 1334. 

For the elastic FEA computation, a prescribed displacement $`\bar{\underline{u}}_{a}^{\#} = [u_x,0,0]`$ is applied on both the ends of the specimen (shown in red) in opposite x-directions.

The von Mises stress coming from this elastic computation $`\bar{\sigma}_{\textrm{VM}}^{\#} = \sigma_y =`$ 170 MPa in the gauge section away from pores. This elastic FEA computation corresponds to $`f(t)`$ = 1. The results can be examined in the paraview file **sigvme_f1.vtu**. These results are also written into the text file **sigvme_f1.txt** for easy input to the plastic correction algorithm.

The load function $`f(t)`$ is chosen to oscillate between +0.47 and -0.47, which leads to the chosen scaled loading $`f(t)\bar{\sigma}_{\textrm{VM}}^{\#}`$ = 80 MPa in the gauge section.

<img src="https://github.com/user-attachments/assets/6817691d-6c80-4f03-90f0-4675df38fa8b" width="792px" height="417px">

Near the pores, the chosen scaled loading $`f(t)\bar{\sigma}_{\textrm{VM}}^{\#}`$ goes up to 350 MPa, which represents a stress concentration factor of $`\sim`$4.4. This is shown in the following image:

<img src="https://github.com/user-attachments/assets/ddecd1aa-3d68-4926-a722-5ff0db5dc2cc" width="423.5px" height="339.5px">

The full-field $`\bar{\sigma}_{\textrm{VM}}^{\#}`$ is input as an array (given by the variable _sig_vm_e_ in **run_model_fullmesh.py**) to the plastic corrector algorithm. The cumulative plastic strain $`p`$ is obtained as output, for all the time-steps defined in the load function $`f(t)`$.


## Accuracy of the plastic correction algorithm

We show here the accuracy of the plastic corrector for $`\Delta p = p^{\textrm{max}}_{\textrm{cycle}} - p^{\textrm{min}}_{\textrm{cycle}}`$ in the 20<sup>th</sup> cycle, by comparing it to a full elasto-plastic FEA computation, which serves as the reference. Users are invited to refer to the full publication [[1]](#1) for more details on the plastic correction algorithm and error analysis on a wider range of BVPs.

<img src="https://github.com/user-attachments/assets/1fba19a5-4a06-4d8a-9d92-8fa05f8a5d88" width="752.6px" height="794.6px">

## References

<a id="1">[1]</a> 
Palchoudhary, Abhishek, Simone Peter, Vincent Maurel, Cristian Ovalle, and Pierre Kerfriden. 
A plastic correction algorithm for full-field elasto-plastic finite element simulations: 
critical assessment of predictive capabilities and improvement by machine learning. 
arXiv preprint arXiv:2402.06313 (2024).

<a id="2">[2]</a> 
M. S. Alnaes, J. Blechta, J. Hake, A. Johansson, B. Kehlet, A. Logg, C. Richardson, J. Ring, M. E. Rognes and G. N. Wells. 
The FEniCS Project Version 1.5, Archive of Numerical Software 3 (2015). 
[doi.org/10.11588/ans.2015.100.20553]

<a id="3">[3]</a>
Le, V.-D., Saintier, N., Morel, F., Bellett, D., Osmond, P. Investigation of the
effect of porosity on the high cycle fatigue behaviour of cast al-si alloy by x-
ray micro-tomography. International Journal of Fatigue 106, 24–37 (2018) 
https://doi.org/10.1016/j.ijfatigue.2017.09.012
