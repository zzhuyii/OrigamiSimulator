%%%%% Sequentially Working Origami Multi-Physics Simulator (SWOMPS)  %%%%%%
%
% Authors: Yi Zhu, and Evgueni T. Filipov
%
% Discription: This code package implement a bar and hinge model based 
% simulator for active origami structures with multi-physics based 
% actuation mechanisms. The code package can capture both the mechanical 
% behavior and the heat transfer aspect. The implementation is versatile
% and has the following features:
%
% (1) Provides 5 different loading solvers of active origami. They are: 
%     Newton-Raphson method, displacement controlled method, modified 
%     generazlied displacement controlled method, self-stress folding, and
%     thermal folding method.
% (2) Allows users to create arbitrary number and sequence of the five
%     loading methods. Users can stop the solver at specified increments
%     and switch between different solvers or edit origami systems during 
%     within the increment easily.
% (3) Simulate electro-thermo-mechanically coupled actuation of origami.
% (4) Simulate inter-panel contact of origami systems.
% (5) Simulate the compliant creases explicitly with novel bar and hinge
%     model meshing schemes.
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% origami simulator, the origami contact, compliant crease, electro-thermal
% model presented in this package. 
%
% Reference:
% [1] Y. Zhu, E. T. Filipov (2021). 'Sequentially Working Origami Multi-
%     Physics Simulator (SWOMPS): A Versatile Implementation' (submitted)
% [2] Y. Zhu, E. T. Filipov (2021). 'Rapid Multi-Physic Simulation for 
%     Electro-Thermal Origami Robotic Systems' (submitted)
% [3] Y. Zhu, E. T. Filipov (2020). 'A Bar and Hinge Model for Simulating 
%     Bistability in Origami Structures with Compliant Creases' Journal of 
%     Mechanisms and Robotics, 021110-1. 
% [4] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' Proc. R. Soc. A, 475: 
%     20190366.       
% [5] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [6] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [7] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' Proc. R. Soc. A 473: 
%     20170348. 
% [8] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%% Sequentially Working Origami Multi-Physics Simulator (SWOMPS)  %%%%%%

%% Initialize the solver
clear all;
clc;
close all;


%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing. First, we define the parameters that will be used
% in the modeling.

% thickness of two layers of material 
t1=0.2*10^-6;
t2=0.80*10^-6;

% thickness of panel
tpanel=21*10^-6;

% density of panel
rho=1200;

% residual folding developed after fabrication (degree)
residualFold=10;

% width of crease
W=400*10^-6;

% length of panel
Lpanel=1*10^(-3);

% power input (mW)
qload=10;
qload = 20;
% to fully close the single crease, we can input 20mW power

% create an instance of the solver class
ori=OrigamiSolver;

% define the Geometry
ori.node0=[0 0 0;
      Lpanel+W/2 0 0;
      Lpanel+W/2 0 Lpanel;
      0 Lpanel 0;
      Lpanel+W/2 Lpanel 0;
      Lpanel+W/2 Lpanel  Lpanel;];
  
ori.panel0{1}=[1 2 5 4];
ori.panel0{2}=[2 3 6 5];


% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();

% Plot the results for inspection
ori.displayRange=3*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;

% Compute the meshed geometry
ori.flag2D3D=3;
ori.compliantCreaseOpen=1;
ori.Mesh_Mesh()

% Plot the meshed origami for inspection;
ori.Plot_MeshedOrigami(); 


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelW=W;
ori.diagonalRate=1000; % crease torsion stiffness/ bending stiffness

ori.panelThickVec=[tpanel;tpanel]; 
ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(t1+t2);


%% setup panel contact information

ori.contactOpen=1;
ori.ke=0.0001;
ori.d0edge=40*(10^(-6));
ori.d0center=40*(10^(-6));


%% Assign Thermal Properties

ori.panelThermalConductVec = [1.3;1.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=1500*10^(-6);



%% Setup the loading controller

% applying zero force
nr=ControllerNRLoading;

nr.increStep=10;
nr.tol=10^-6;
nr.iterMax=50;

nr.supp=[1,1,1,1;
      2,1,1,1;
      7,1,1,1;];

loadMag=0;
nr.load=[5,0,0,-loadMag;];
nr.videoOpen=0;

% applying the residual stress
selfFold=ControllerSelfFolding;

% Assign zero strain position for creases during self-folding
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge. This should fold back flat
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);
selfFold.targetRotZeroStrain(3)=pi;

if ori.compliantCreaseOpen==1
    selfFold.supp=[1,1,1,1;
          2,1,1,1;
          12,1,1,1;];
else
    selfFold.supp=[1,1,1,1;
          2,1,1,1;
          7,1,1,1;];
end
   
selfFold.increStep=80;
selfFold.tol=1*10^-9;
selfFold.iterMax=25;
selfFold.videoOpen=0;


ori.loadingController{1}={"NR",nr};
ori.loadingController{2}={"SelfFold",selfFold};


%% Solving the model
ori.Solver_Solve();
