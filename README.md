# OrigamiSimulator

![alt text](https://github.com/zzhuyii/OrigamiSimulator/blob/master/10_Documents_Figures/SWOMPS.gif)

Realistic simulation of active origami structures with multi-physics behaviors. 
The package can capture compliant creases, inter-panel contact, heat transfer, large
folding deformation, and electro-thermally actuated creases. The package provides a wide variety of
loading and simulation methods and allows users to created customizable loading schemes with 
arbitrary number of loading sequences. I am actively adding new capabilities to the package over time. 

## Main Features of the package:

* Loading & simulation methods: 
    * Newton-Raphson method
    * Displacement controlled method
    * Generalized displacement controlled method for applying external forces
    * Changing stress-free folding angle for self-folding
    * Changing environmetnal temperature
    * Applying electro-thermal actuation for active self-folding
    * Constant Average Acceleration Method for dynamic loading
    * Linear Multistep Method for transient heat transfer coupled with origami folding
    
* Allows users to create customizable loading schemes with arbitrary number and 
    sequence of the provided loading methods. 
    
* Simulates compliant creases in active origami and provides automated Meshing 
    code for compliant creases. (**Figure 1**)

* Simulates inter-panel contact induced mechanical behaviors within origami. (**Figure 2**)

* Simulates heat transfer in origami systems and captures the electro-
    thermo-mechanically coupled actuation of active origami creases. (**Figure 3**)

* Simulates the kinematic folding and mechanical load-carrying behaviors of thick origami structures. 
    
<p align="center">
<img src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/10_Documents_Figures/RealisticCrease.png" width="45%" >
</p>

**Figure 1.** The package allows users to simulate origami systems with compliant creases, which are creases with non-negligible width. 
Using the compliant crease meshing provides more realistic geometry and allows the simulator to capture advanced mechanical behaviors such as bistability
 and multi-physics actuation. 

<p align="center">
<img align="center" src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/10_Documents_Figures/InterLockingDevice.PNG" width="60%" >
</p>

**Figure 2.** The package allows users to simulate global inter-panel contact within the origami. This panel contact model is a physics-based 
frictionless contact model and can give the correponding forcese between the contacting panels. 
 
<p align="center">
<img align="center" src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/10_Documents_Figures/GraphAbstract_webpage.png" width="50%" >
</p>
    
**Figure 3.** The package allows users to simulate multi-physics based electro-thermal actuation in the origami creases. The package can capture 
the heat transfer of active origami and can calculate the active crease folding due to changing temperature. 


## Comments on Accuracy in Stiffness Prediction:

This code is based on bar and hinge models for origami structures. The bar and hinge model is a reduced order simulation method
for origami structures. There are different methods to derive the stiffness parameters of a bar and hinge model, and they produce different results. 
In this simulation package, the stiffness parameters are derived analytically through matching the stiffness of bar and hinge models
to theoretical plates (for shear and tension) and pseudo-rigid-body models (for large folding). 
This analytical method does not capture large panel deformations, especially when studying non-rigid-foldable origami patterns. 
For example, most paper-based non-rigid-foldable origami prototypes can have panels with initial curvature and bending, which cannot be captured using bar and hinge models. 
Therefore, this simulation package can potentially over-estimate the stiffness of non-rigid foldable origami systems.
To resolve this, many research use curve-fitting to find the stiffness parameters of bars and rotational springs (bar areas and rotational spring stiffness). 
In this case, curve-fitting is like using softened secant stiffness (while analytically derivation is like using the initial tangent stiffness), 
so that the model can better approximate the behaviors of non-rigid-foldable origami prototypes. 
Users are suggested to try curve-fit the stiffness parameters for better accuracy when studying non-rigid-foldable origami systems. 
However, the simulation results from analytical derivation still provides a fast alternative for finding the trends and understanding the behaviors of origami systems. 

## Efficiency Update (2022-07-11 & 2022-07-19)

Part of the code is vectorized now. I cannot believe that I did not do it previously.
This probably shows the beauty of bar and hinge models; you can follow bad coding habbit and still get reasonalbly fast results. 
Thankfully, one of my wafer is stuck in the machine and the lab staff is crazily busy recently so he cannot remove it for me. 
This prevents me from doing micro-fab so I have an entire weekend plus a couple week days to sit down and vectorize my previous codes. 
In addition to vectorizing the code for mechanical simulation, part of the code for thermal simulaiton is streamlined to reduce redundant calculation. 

<p align="center">
<img align="center" src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/10_Documents_Figures/Efficiency Update.jpg" width="100%" >
</p>
    
**Figure 4.** 
For the "Example06_FlowerSelfFold.m" the new version is roughly five times faster than the previous code. 
For the simulation of the SWOMPS logo, the new code is about three times faster than the previous code. 


## Using the Code:

PLEASE ADD THE "00_SourceCode" IN TO THE PATH. 


## Acknowledgement: 

I would like to acknowledge the prior works from Ke Liu and Glaucio H. Paulino 
for non-rigid origami simulators.  Their works pave the ground for the development
of this package.


## Reference:

1. Y. Zhu, E. T. Filipov (**2024**) Large-scale modular and uniformly thick origami-inspired
    adaptable and load-carrying structures, Nature Communications, 15, 2353.

1. Y. Zhu, E. T. Filipov (**2021**) Sequentially Working Origami Multi-Physics Simulator 
    (SWOMPS): A Versatile Implementation, IDETC-CIE 2021, DETC2021-68042.

1. Y. Zhu, E. T. Filipov (**2021**) Rapid Multi-Physics Simulation for Electro-Thermal 
    Origami Systems.  *International Journal of Mechanical Sciences*, 202-203, 106537.     

1. Y. Zhu, E. T. Filipov (**2020**). A Bar and Hinge Model for Simulating Bistability 
    in Origami Structures With Compliant Creases, *Journal of Mechanisms and Robotics*, 
    12, 021110-1.
    
1.	Y. Zhu, E. T. Filipov (**2019**). An Efficient Numerical Approach for Simulating 
    Contact in Origami Assemblage, *Proceedings of the Royal Society A*, 475, 20190366.
    
1.	Y. Zhu, E. T. Filipov (**2019**). Simulating compliant crease origami with a bar and
    hinge model. *IDETC/CIE 2019*, 97119. 

1.	K. Liu, G. H. Paulino (**2018**). Highly efficient nonlinear structural analysis of 
    origami assemblages using the MERLIN2 software, *Origami^7*.
    
1.	K. Liu, G. H. Paulino (**2017**). Nonlinear mechanics of non-rigid origami – An efficient
    computational approach, *Proceedings of the Royal Society A*, 473, 20170348.
    
1.	K. Liu, G. H. Paulino (**2016**). MERLIN: A MATLAB implementation to capture highly 
    nonlinear behavior of non-rigid origami, *Proceedings of IASS Annual Symposium 2016*. 
