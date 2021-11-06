# OrigamiSimulator

![alt text](https://github.com/zzhuyii/OrigamiSimulator/blob/master/04_Documents_Figures/SWOMPS.gif)

Realistic simulation of active origami structures with multi-physics behaviors. 
The package can capture compliant creases, inter-panel contact, heat transfer, large
folding deformation, and electro-thermally actuated creases. The package provides five
different loading methods adn allows users to created customizable loading schemes with 
arbitrary number and sequence of the provided loading methods. 

Please check the package on our group website or my personal website:

https://drsl.engin.umich.edu/software/swomps-package/

https://sites.google.com/view/yi-zhu/research/origami-simulator/

## Main Features of the package:

* Provides five different loading simulation methods, inlcuding Newton-Raphson method,
    Displacement controlled method, Generalized displacement controlled method,
    self folding method, and electro-thermal actuation method. 
    
* Allows users to create customizable loading schemes with arbitrary number and 
    sequence of the provided five loading methods. 
    
* Simulates compliant creases in active origami and provides automated Meshing 
    code for compliant creases. (**Figure 1**)

* Simulates inter-panel contact induced mechanical behaviors within origami. (**Figure 2**)

* Simulates heat transfer in origami systems and captures the electro-
    thermo-mechanically coupled actuation of active origami creases. (**Figure 3**)
    
<p align="center">
<img src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/04_Documents_Figures/RealisticCrease.png" width="45%" >
</p>

**Figure 1.** The package allows users to simulate origami systems with compliant creases, which are creases with non-negligible width. 
Using the compliant crease meshing provides more realistic geometry and allows the simulator to capture advanced mechanical behaviors such as bistability
 and multi-physics actuation. 

<p align="center">
<img align="center" src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/04_Documents_Figures/InterLockingDevice.PNG" width="60%" >
</p>

**Figure 2.** The package allows users to simulate global inter-panel contact within the origami. This panel contact model is a physics-based 
frictionless contact model and can give the correponding forcese between the contacting panels. 
 
<p align="center">
<img align="center" src="https://github.com/zzhuyii/OrigamiSimulator/blob/master/04_Documents_Figures/GraphAbstract_webpage.png" width="50%" >
</p>
    
**Figure 3.** The package allows users to simulate multi-physics based electro-thermal actuation in the origami creases. The package can capture 
the heat transfer of active origami and can calculate the active crease folding due to changing temperature. 

## Using the Code:

PLEASE ADD THE "00_SourceCode" IN TO THE PATH. For standard mechanical simulation 
of origami, please check the selected simulaiton example from the JMR paper and
PRSA paper in "01_MechanicalLoadingExample". For simulation of the folding electro-thermally 
active origami, please check the example codes in "02_ThermalLoadingExample" associated with the IJMS paper. Additional 
tutuorial examples are presented in the "03_SampleCode_Tutorial_IDETC".


## Acknowledgement: 

We would like to acknowledge the prior works from Ke Liu and Glaucio H. Paulino 
for non-rigid origami simulators.  Their works pave the ground for the development
of this package.


## Reference:


1. Y. Zhu, E. T. Filipov (**2021**) SEQUENTIALLY WORKING ORIGAMI MULTI-PHYSICS SIMULATOR 
    (SWOMPS): AVERSATILE IMPLEMENTATION (Accepted to IDETC 2021)

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
