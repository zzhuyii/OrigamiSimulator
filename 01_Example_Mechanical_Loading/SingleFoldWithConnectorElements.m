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
clear;clc;close all;
ori=OrigamiSolver;

%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% thickness of panel
tpanel=21*10^-6;

% thickness of crease
tcrease=2*10^-6;

% length of panel
Lpanel=1*10^(-3);

% Width of crease
W=0.1*10^(-3);

ori.node0=[0 0 0;
      Lpanel 0 0;
      Lpanel Lpanel 0;
      0 Lpanel 0;
      Lpanel 0 0;
      2*Lpanel 0 0;
      2*Lpanel Lpanel 0;
      Lpanel Lpanel 0;
      ];
  
ori.panel0{1}=[1 2 3 4];
ori.panel0{2}=[5 6 7 8];


% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=3*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;



%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;

ori.compliantCreaseOpen=0;

% Compute the meshed geometry
ori.Mesh_Mesh()
% Plot the meshed origami for inspection;
ori.Plot_MeshedOrigami(); 


%% Assign Mechanical Properties & Initialize

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[500*10^(-6);tpanel;500*10^(-6)]; 
ori.panelW=W;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=1000;
ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(tcrease);


%% Add Connector
ori.connectorNode(1,:)=[2,5];
ori.connectorNode(2,:)=[3,8];

ori.connectorK=1000000;
ori.connectorOpen=1;

%% Add Hinges
crossBarArea=2*10^-8;
hingeStiff=0.000001;

% State the additional hinges
ori.barType(17)=1;
ori.barConnect(17,:)=[5,8];
ori.barArea(17)=crossBarArea;
ori.barLength(17)=Lpanel;        
ori.sprK(17)=hingeStiff;
ori.sprIJKL(17,:)=[9,5,8,10];


% This code initialize the solver
ori.Solver_Solve();
ori.continuingLoading=1;

% The initialization will change some of the proeprties
% Need to reasign the properties
ori.barType(17)=1;
ori.barConnect(17,:)=[5,8];
ori.barArea(17)=crossBarArea;
ori.barLength(17)=Lpanel;        
ori.sprK(17)=hingeStiff;
ori.sprIJKL(17,:)=[9,5,8,10];



%% Loading steps
% applying the gravity loading
nr=ControllerNRLoading;

nr.increStep=50;
nr.tol=10^-6;
nr.iterMax=50;

nr.supp=[1,1,1,1;
      9,1,1,1;
      4,1,1,1;];

loadMag=1*10^-4;
nr.load=[6,0,0,loadMag;
      7,0,0,loadMag;];
  
nr.videoOpen=1;
ori.loadingController{1}={"NR",nr};
ori.Solver_Solve();

