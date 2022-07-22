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

tic
%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% size of the box
a=10*10^-3;
b=1.5*10^-3;
l=30*10^-3;


% thickness of two layer
t1=0.02*10^-3;
t2=0.02*10^-3;

% thickness of panel
tpanel=1*10^-3;

% width of crease
W=2*10^-3;

% power input (W)
qload=0.2;

ori.node0=[0 -a 0;
           0 0 0;
           0 a 0;
           2*l+b -a 0;
           2*l-b 0 0;
           2*l+b a 0;
           3*l+b -a 0;
           3*l-b 0 0;
           3*l+b a 0;
           4*l+b -a 0;
           4*l-b 0 0;
           4*l+b a 0;
           6*l -a 0;
           6*l 0 0;
           6*l a 0;
           ];
  
ori.panel0{1}=[1 4 5 2];
ori.panel0{2}=[2 5 6 3];
ori.panel0{3}=[4 7 8 5];
ori.panel0{4}=[5 8 9 6];
ori.panel0{5}=[7 10 11 8];
ori.panel0{6}=[8 11 12 9];
ori.panel0{7}=[10 13 14 11];
ori.panel0{8}=[11 14 15 12];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();

% Plot the results for inspection
ori.viewAngle1=10;
ori.viewAngle2=15;
ori.displayRange=150*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis
ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;

%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;
ori.creaseWidthVec(4)=W;
ori.creaseWidthVec(6)=W;
ori.creaseWidthVec(10)=W;
ori.creaseWidthVec(9)=W;
ori.creaseWidthVec(11)=W;
ori.creaseWidthVec(15)=W;
ori.creaseWidthVec(14)=W;
ori.creaseWidthVec(16)=W;
ori.creaseWidthVec(20)=W;


% Compute the meshed geometry
ori.Mesh_Mesh()




ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

ori.panelThickVec=tpanel*ones(16,1); 
ori.panelW=W;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=1000;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);

ori.creaseThickVec(3)=t1+t2;
ori.creaseThickVec(4)=t1+t2;
ori.creaseThickVec(6)=t1+t2;
ori.creaseThickVec(10)=t1+t2;
ori.creaseThickVec(9)=t1+t2;
ori.creaseThickVec(11)=t1+t2;
ori.creaseThickVec(15)=t1+t2;
ori.creaseThickVec(14)=t1+t2;
ori.creaseThickVec(16)=t1+t2;
ori.creaseThickVec(20)=t1+t2;


%% Assign Thermal Properties

ori.panelThermalConductVec = 0.3*ones(16,1); 

ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=0.02; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=250;
thermal.tol=5*10^-5; 

thermal.supp=[1,1,1,1;
      26,0,1,1;
      8,1,0,1;
      31,0,0,1;];

thermal.thermalBoundaryPanelVec=[16];
thermal.roomTempNode=[];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
deltaAlpha=1000*10^(-6);
thermal.deltaAlpha(3)= -deltaAlpha;
thermal.deltaAlpha(4)= 0;
thermal.deltaAlpha(6)= -deltaAlpha;
thermal.deltaAlpha(10)= 0;
thermal.deltaAlpha(9)= deltaAlpha;
thermal.deltaAlpha(11)= deltaAlpha;
thermal.deltaAlpha(15)= 0;
thermal.deltaAlpha(14)= -deltaAlpha;
thermal.deltaAlpha(16)= -deltaAlpha;
thermal.deltaAlpha(20)= 0;


thermal.Emat1=2*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=t1;
thermal.tmat2=t2;
thermal.videoOpen=0;

% the target loading of crease heating
thermal.targetCreaseHeating=[
    3,qload;
    4,qload;
    6,qload;
    10,qload;
    9,qload;
    11,qload;
    15,qload;
    14,qload;    
    16,qload;
    20,qload;];

ori.loadingController{1}={"ElectroThermal",thermal};
ori.Solver_Solve();

toc
