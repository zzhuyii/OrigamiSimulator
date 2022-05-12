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

W=0.4;
L=1;
L1=1;
L2=1;
tc=0.001;
tp1=0.01;
tp2=0.01;

ori.node0=[0 0 0;
      L1+W/2 0 0;
      L2+L1+W 0 0;
      0 L 0;
      L1+W/2 L 0;
      L2+L1+W L 0;];
  
ori.panel0{1}=[1 2 5 4];
ori.panel0{2}=[2 3 6 5];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=3; % plotting range
ori.displayRangeRatio=0.1; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[tp1;tp2];
ori.panelW=W;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=tc;

%% setup panel contact information

ori.contactOpen=0;
ori.ke=0.0001;
ori.d0edge=40*(10^(-6));
ori.d0center=40*(10^(-6));

%% Assign Thermal Properties

ori.panelThermalConductVec = [1;1]; 
ori.creaseThermalConduct=1;
ori.envThermalConduct=0;

% thickness of the submerged environment at RT
ori.t2RT=1;

ori.envLayer=1;
ori.thermalDissipation=20/180*3.14;
ori.RT=0;



%% Thermal Loading

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=1;
thermal.tol=5*10^-7; 
thermal.iterMax=10;

thermal.supp=[1,1,1,1;
      2,1,1,1;
      3,1,1,1;
      4,1,1,1;
      9,1,1,1;
      10,1,1,1;
      11,1,1,1;
      12,1,1,1;];

thermal.thermalBoundaryPanelVec=[-1];
thermal.roomTempNode=[1;4;6;7];

thermal.deltaAlpha(3)=0; 
thermal.Emat1=1;
thermal.Emat2=1;
thermal.tmat1=1;
thermal.tmat2=1;
thermal.videoOpen=0;

% the target loading of crease heating
thermal.targetCreaseHeating=[3,10*W*L*tc];

ori.loadingController{1}={"ThermalLoading",thermal};


%% Solving the model
ori.Solver_Solve();

T=thermal.temperatureHis;
Tave=1/3*(T(9)+T(10)+T(11))
Tmax=T(11)

