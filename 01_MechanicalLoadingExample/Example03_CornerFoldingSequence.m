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
clear all;clc;close all;
ori=OrigamiSolver;

%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% Define the nodal coordinate before meshing
a=20*10^(-3);
ori.node0=[0 0 0;
      2*a 0 0;
      3*a 0 0;
      0 2*a 0; %4
      2*a 2*a 0;
      3*a 2*a 0;
      4*a 2*a 0; %7
      0 3*a 0;
      2*a 3*a 0;
      3*a 3*a 0;
      4*a 3*a 0;];
  
ori.panel0{1}=[1 2 5 4];
ori.panel0{2}=[2 3 6 5];
ori.panel0{3}=[4 5 9 8];
ori.panel0{4}=[5 6 10];
ori.panel0{5}=[6 7 11 10];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=4*10^(-3);
ori.creaseWidthVec(4)=4*10^(-3);
ori.creaseWidthVec(7)=4*10^(-3);
ori.creaseWidthVec(12)=4*10^(-3);
ori.creaseWidthVec(13)=4*10^(-3);

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=100*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis


% We need to manually adjust some nodal coordinates for this case
creaseW=4*10^(-3);
ori.newNode(10,1)=ori.newNode(10,1)-0.5*creaseW;
ori.newNode(11,1)=ori.newNode(11,1)-0.5*creaseW;
ori.newNode(23,1)=ori.newNode(23,1)-0.25*creaseW;

ori.newNode(14,1)=ori.newNode(14,1)+0.5*creaseW;
ori.newNode(15,1)=ori.newNode(15,1)+0.5*creaseW;
ori.newNode(26,1)=ori.newNode(26,1)+0.25*creaseW;

ori.newNode(29,1)=ori.newNode(29,1)+0.5*creaseW;
ori.newNode(30,1)=ori.newNode(30,1)+0.5*creaseW;
ori.newNode(31,1)=ori.newNode(31,1)+0.5*creaseW;
ori.newNode(36,1)=ori.newNode(36,1)+0.25*creaseW;

ori.newNode(19,1)=ori.newNode(19,1)+0.5*creaseW;
ori.newNode(16,1)=ori.newNode(16,1)+0.5*creaseW;

ori.newNode(19,2)=ori.newNode(19,2)-0.5*creaseW;
ori.newNode(16,2)=ori.newNode(16,2)+0.5*creaseW;

ori.newNode(18,2)=ori.newNode(18,2)-0.5*creaseW;
ori.newNode(17,2)=ori.newNode(17,2)+0.5*creaseW;

ori.newNode(29,2)=ori.newNode(29,2)+0.25*creaseW;
ori.newNode(30,2)=ori.newNode(30,2)-0.25*creaseW;

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[1;1;1;1;1]*200*10^(-6); 
ori.panelW=4*10^-3;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=100*10^(-6);
ori.creaseThickVec(4)=100*10^(-6);
ori.creaseThickVec(7)=100*10^(-6);
ori.creaseThickVec(12)=100*10^(-6);
ori.creaseThickVec(13)=100*10^(-6);

%% setup panel contact information
ori.contactOpen=1;
ori.ke=0.02;
ori.d0edge=1.5*10^-3;
ori.d0center=1.5*10^-3;

%% Setup the loading controller

% define the self folding step
selfFold=ControllerSelfFolding;

% Assign zero strain position for creases during self-folding
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);

selfFold.targetRotZeroStrain(3)=pi+0.5*pi;
selfFold.targetRotZeroStrain(4)=pi+0.5*pi;

selfFold.supp=[1,1,1,1;
          2,1,1,1;
          3,1,1,1;
          4,1,1,1;];
 
selfFold.increStep=40;
selfFold.tol=10^-5;
selfFold.iterMax=50;


%% Multi step loading
% Here we demonstrate another way to achieve the multi-step loading using
% the "contiuingLoading" option in the package.

% solve the system
ori.loadingController{1}={"SelfFold",selfFold};
ori.Solver_Solve();

% A second step of loading
ori.continuingLoading=1;
selfFold.targetRotZeroStrain(7)=pi+0.5*pi;
ori.Solver_Solve();

% A third step of loading
selfFold.targetRotZeroStrain(12)=pi+pi;
ori.Solver_Solve();
