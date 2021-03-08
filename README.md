# OrigamiSimulator
Simulating active origami structures with multi-pyhsics behaviors including 
compliant creases, inter-panel contact, and electro-thermally actuated creases. 


Main Features of the package:

(1) Provides five different loading simulation methods: Newton-Raphson method,
    Displacement controlled method, Generalized displacement controlled method,
    self folding method, and electro-thermal actuation method. 
    
(2) Allows users to create customized loading schemes with arbitrary number and 
    sequence of the providev five loading methods. 
    
(3)	Simulates compliant creases in active origami and provides automated Meshing 
    code for compliant creases.

(4)	Simulates inter-panel contact induced mechanical behaviors within origami.

(5) Simulates heat transfer in origami systems, and can capture the electro-
    thermo-mechanically coupled actuation of active origami creases. 


Using the Code:

PLEASE ADD THE "00_SourceCode" IN TO THE PATH. For standard mechanical simulation 
of origami, please check the selected simulaiton example from the JMR paper and
PRSA paper in "01_MechanicalLoadingExample".

For simulation of the folding electro-thermally active origami, please check the
example codes in "02_ThermalLoadingExample".


Acknowledgement: 

We would like to acknowledge the prior works from Ke Liu and Glaucio H. Paulino 
for non-rigid origami simulators.  Their works pave the ground for the development
of this package.


Reference:
[1] Y. Zhu, E. T. Filipov (2021) SEQUENTIALLY WORKING ORIGAMI MULTI-PHYSICS SIMULATOR 
    (SWOMPS): AVERSATILE IMPLEMENTATION (submitted)

[2] Y. Zhu, E. T. Filipov (2021) Rapid Multi-Physics Simulation for Electro-Thermal 
    Origami Systems (submitted)

[3] Y. Zhu, E. T. Filipov (2020). A Bar and Hinge Model for Simulating Bistability 
    in Origami Structures With Compliant Creases, Journal of Mechanisms and Robotics, 
    12, 021110-1.
    
[4]	Y. Zhu, E. T. Filipov (2019). An Efficient Numerical Approach for Simulating 
    Contact in Origami Assemblage, Proceedings of the Royal Society A, 475, 20190366.
    
[5]	Y. Zhu, E. T. Filipov (2019). Simulating compliant crease origami with a bar and
    hinge model. IDETC/CIE 2019, 97119. 

[6]	K. Liu, G. H. Paulino (2018). Highly efficient nonlinear structural analysis of 
    origami assemblages using the MERLIN2 software, Origami^7.
    
[7]	K. Liu, G. H. Paulino (2017). Nonlinear mechanics of non-rigid origami – An efficient
    computational approach, Proceedings of the Royal Society A, 473, 20170348.
    
[8]	K. Liu, G. H. Paulino (2016). MERLIN: A MATLAB implementation to capture highly 
    nonlinear behavior of non-rigid origami, Proceedings of IASS Annual Symposium 2016. 
