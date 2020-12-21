%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Yi Zhu, and Evegueni T. Filipov
%
% Discription: This code package implement a bar and hinge model based 
% simulator for active origami structures. The code package can capture 
% the following behaviors of active origami systems;
%
% (1) Simulate panel contact in origami;
% (2) Simulate electro-thermal actuation for folding origami;
% (3) Provide 3 different solver for large deformation loading;
% (4) Provide elastic support;
% (5) Provide compliant crease bar and hinge model for origami;
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% origami simulator, the origami contact, compliant crease, electro-thermal
% model presented in this package. 
%
% Reference:
% [1] Y. Zhu, E. T. Filipov (2020). 'Rapid Multi-Physic Simulation for 
%     Electro-Thermal Origami Robotic Systems' (submitted)
% [2] Y. Zhu, E. T. Filipov (2020). 'A Bar and Hinge Model for Simulating 
%     Bistability in Origami Structures with Compliant Creases' Journal of 
%     Mechanisms and Robotics, 021110-1. 
% [3] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' Proc. R. Soc. A, 475: 
%     20190366.       
% [4] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [5] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [6] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' Proc. R. Soc. A 473: 
%     20170348. 
% [7] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize the solver
clear all;;clc;close all;
ori=OrigamiSolver;

%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% Define the nodal coordinate before meshing
a=20*10^(-3);
ori.node0=[0 0 0;
      0 a 0;
      0 a*2.6 0;
      0 a*3.4 0;
      a 0 0;
      a a 0;
      a a*2.6 0;
      a a*3.4 0;];
  
% Define the panel connectivity before meshing
ori.panel0{1}=[5 6 2 1];
ori.panel0{2}=[6 7 3 2];
ori.panel0{3}=[7 8 4 3];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model
% Define the crease width 
ori.creaseWidthMat=zeros(ori.oldCreaseNum,1);
ori.creaseWidthMat(3)=1.5*10^(-3);
ori.creaseWidthMat(6)=1.5*10^(-3);

% Compute the meshed geometry
ori.Mesh_CompliantCreaseGeometry()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=60*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickMat=[1;1;1]*500*10^(-6);
ori.panelW=1.5*10^-3;

ori.creaseThickMat=zeros(ori.oldCreaseNum,1);
ori.creaseThickMat(3)=90*10^(-6);
ori.creaseThickMat(6)=90*10^(-6);

%% setup panel contact information
ori.contactOpen=1;
ori.ke=0.002;
ori.d0edge=1.5*10^-3;
ori.d0center=0.7*1.5*10^-3;

%% Setup the loading controller
% define the self folding step
selfFold=ControllerSelfFolding;

% Assign zero strain position for creases during self-folding
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);

ratio=0.9;
selfFold.targetRotZeroStrain(3)=pi+ratio*pi;
selfFold.targetRotZeroStrain(6)=pi-ratio*pi;

selfFold.supp=[1,1,1,1;
               2,1,0,1;
               3,0,0,1;
               4,0,0,1;];
selfFold.increStep=20;
selfFold.tol=10^-6;
selfFold.iterMax=50;

% define a NR loading step
nr=ControllerNRLoading;

nr.supp=[1,1,1,1;
         2,1,0,1;
         3,0,0,1;
         4,0,0,1;];
loadForce=0.5*10^(-3);
nr.load=[10,0,0,-loadForce;
         11,0,0,-loadForce;];
nr.increStep=50;
nr.tol=10^-6;
nr.iterMax=50;


ori.loadingController{1}={"SelfFold",selfFold};
ori.loadingController{2}={"NR",nr};

%% Solving the model
ori.Solver_Solve();

