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
a=50*10^(-3);
b=50*10^(-3);
z=7.4668*10^(-3);

% Define the nodal coordinate before meshing
ori.node0=[0 0 z;
      -b 0 0;
      0 -b 0;
      b 0 0;
      0 b 0;
      -a -a -z;
      a -a -z;
      a a -z;
      -a a -z;];

% Define the panel connectivity before meshing
ori.panel0{1}=[2 1 5 9];
ori.panel0{2}=[1 4 8 5];
ori.panel0{3}=[6 3 1 2];
ori.panel0{4}=[3 7 4 1];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model
% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(2)=6*10^(-3);
ori.creaseWidthVec(3)=6*10^(-3);
ori.creaseWidthVec(5)=6*10^(-3);
ori.creaseWidthVec(10)=6*10^(-3);

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.displayRange=60*10^(-3); % plotting range
ori.displayRangeRatio=1; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2852*10^6; 
ori.creaseE=25*10^6;
ori.panelPoisson=0.4;
ori.creasePoisson=0.4; 
ori.panelThickVec=[1;1;1;1]*2.2*10^(-3);
ori.panelW=6*10^-3;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(2)=1*10^(-3);
ori.creaseThickVec(3)=1*10^(-3);
ori.creaseThickVec(5)=1*10^(-3);
ori.creaseThickVec(10)=1*10^(-3);


%% Setup the loading controller

% we use the displacement controlled method for the loading
dc=ControllerDCLoading;
dc.supp=[2,1,1,1;
         5,1,1,1;
         11,1,1,1;
         16,1,1,1;];     
loadForce=3;
dc.load=[29,0,0,loadForce;
      30,0,0,loadForce;
      31,0,0,loadForce;
      32,0,0,loadForce;];  
dc.increStep=50; 
dc.tol=10^-6; 
dc.selectedRefDisp=[29,3]; % use node 29 in z direction as reference

% we then use the MGDCM for the unloading
mgdcm=ControllerMGDCMLoading;
mgdcm.supp=[2,1,1,1;
         5,1,1,1;
         11,1,1,1;
         16,1,1,1;];     
loadForce=-1;
mgdcm.load=[29,0,0,loadForce;
      30,0,0,loadForce;
      31,0,0,loadForce;
      32,0,0,loadForce;];  
mgdcm.increStep=10; 
mgdcm.tol=10^-6; 

ori.loadingController{1}={"DC",dc};
ori.loadingController{2}={"MGDCM",mgdcm};


%% Solving the origami
ori.Solver_Solve();

