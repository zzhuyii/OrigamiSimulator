# OrigamiSimulator
Simulating the mechanical behavior of origami structures with contact, compliant creases, and multi-physics 
behaviros associated with the elctro-thermal actuation 

PLEASE ADD THE "00_SourceCode" IN TO THE PATH.


Main Features of the Code:



(1)	Simulate Compliant Creases in Origami Patterns

(2)	Automated Meshing for Compliant Creases

(3)	Simulate Contact in Origami patterns

(4)	Provide 4 solvers for mechanical loading of origami including:
    Newton-Raphson loading method, Displacement controlled method, Generalized displacement controlled method,
    and a Newton-Raphson solver for self folding

(5)	Non-rigid support available

(6) Simulate Heat Transfer in Origami Systems

(7) Provide a solver for thermal loading: Solve the equilibrium configuration under applied heating power



Using the Code:

For standard mechanical simulation of origami, please check the selected simulaiton example
from the JMR paper and PRSA paper in "01_MechanicalLoadingExample".

For simulation of the folding electro-thermally active origami, please check the
example codes in "02_ThermalLoadingExample".


Acknowledgement: 

We would like to acknowledge the prior works from Ke Liu and Glaucio H. Paulino 
for non-rigid origami simulators.  Their works pave the ground for the development of this package.


Reference:

[1] Y. Zhu, E. T. Filipov (2020) Rapid Multi-Physics Simulation for Electro-Thermal Origami Robotic Systems (submitted)

[2] Y. Zhu, E. T. Filipov (2020). A Bar and Hinge Model for Simulating Bistability in Origami Structures With Compliant Creases,
    Journal of Mechanisms and Robotics, 12, 021110-1.
    
[3]	Y. Zhu, E. T. Filipov (2019). An Efficient Numerical Approach for Simulating Contact in Origami Assemblage,
    Proceedings of the Royal Society A, 475, 20190366.
    
[4]	Y. Zhu, E. T. Filipov (2019). Simulating compliant crease origami with a bar and hinge model. IDETC/CIE 2019, 97119. 

[5]	K. Liu, G. H. Paulino (2018). Highly efficient nonlinear structural analysis of origami assemblages using the MERLIN2 software, 
    Origami^7.
    
[6]	K. Liu, G. H. Paulino (2017). Nonlinear mechanics of non-rigid origami – An efficient computational approach,
    Proceedings of the Royal Society A, 473, 20170348.
    
[7]	K. Liu, G. H. Paulino (2016). MERLIN: A MATLAB implementation to capture highly nonlinear behavior of non-rigid origami, 
    Proceedings of IASS Annual Symposium 2016. 
