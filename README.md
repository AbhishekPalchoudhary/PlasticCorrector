This algorithm pertains to the numerical implementation of a new local plastic correction algorithm that is aimed to
accelerate elasto-plastic finite element (FE) simulations for structural problems
exhibiting localised plasticity (around e.g. notches, geometrical defects). The pro-
posed method belongs to the category of generalised multi-axial Neuber-type
methods, which process the results of an elastic prediction point-wise in order to
calculate an approximation of the elasto-plastic solution. The publication associated to this algorithm can be found here [[2]](#2).

We demonstrate the usage of the plastic corrector with the following example boundary value problem (BVP).

A specimen with a sub-volume of pores explicitly meshed is considered. These pores arise due to a casting manufacturing process, and information on their geometry was obtained from computed tomography.

The material parameters (identified by other authors [[1]](#1)) correspond to an AlSi7Mg0.3 alloy:

Elasticity parameters:

$`E =`$ 75500 MPa, 
$`\sigma_y =`$ 170 MPa. 

Isotropic and kinematic hardening parameters:

$`b =`$ 19, 
$`Q =`$ 20 MPa, 
$`C =`$ 127499 MPa, 
$`D =`$ 1334. 

For the elastic FEA computation, a prescribed displacement $`\bar{\underline{u}}_{a}^{\#} = [u_x,0,0]`$ is applied on both the ends of the specimen (shown in red) in opposite x-directions.

The von Mises stress coming from this elastic computation $`\bar{\sigma}_{\textrm{VM}}^{\#} = \sigma_y =`$ 170 MPa in the gauge section away from pores.

$`f(t)`$ is chosen to oscillate between +0.47 and -0.47, which leads to $`f(t)\bar{\sigma}_{\textrm{VM}}^{\#}`$ = 80 MPa in the gauge section.

<img src="https://github.com/user-attachments/assets/6817691d-6c80-4f03-90f0-4675df38fa8b" width="792px" height="417px">

Near the pores, $`f(t)\bar{\sigma}_{\textrm{VM}}^{\#}`$ goes up to 350 MPa, which represents a stress concentration factor of $`\sim`$4.4.

<img src="https://github.com/user-attachments/assets/ddecd1aa-3d68-4926-a722-5ff0db5dc2cc" width="423.5px" height="339.5px">

The full-field $`\bar{\sigma}_{\textrm{VM}}^{\#}`$ is input as an array (given by the variable _sig_vm_e_ in **run_model.py**) to the plastic corrector algorithm. Elasto-plastic variables like the cumulative plastic strain $`p`$ is obtained as output.

We show here the accuracy of the plastic corrector for $\Delta p$ in the 20<sup>th</sup> cycle, by comparing it to a full elasto-plastic FEA computation, which serves as the reference. Users are invited to refer to the full publication [[2]](#2) for more details on the plastic correction algorithm and error analysis on a wider range of BVPs.

<img src="https://github.com/user-attachments/assets/1fba19a5-4a06-4d8a-9d92-8fa05f8a5d88" width="752.6px" height="794.6px">

## References
<a id="1">[1]</a> 
Le, V.-D., Saintier, N., Morel, F., Bellett, D., Osmond, P.: Investigation of the
effect of porosity on the high cycle fatigue behaviour of cast al-si alloy by x-
ray micro-tomography. International Journal of Fatigue 106, 24–37 (2018) https:
//doi.org/10.1016/j.ijfatigue.2017.09.012

<a id="2">[2]</a>
Palchoudhary, Abhishek, Simone Peter, Vincent Maurel, Cristian Ovalle, and Pierre Kerfriden. 
"A plastic correction algorithm for full-field elasto-plastic finite element simulations: 
critical assessment of predictive capabilities and improvement by machine learning." 
arXiv preprint arXiv:2402.06313 (2024).
