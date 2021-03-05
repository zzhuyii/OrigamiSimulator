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

%% Input data for Defining the Geometry
% This section of code is used to generate the geometry of structure
% Input files are still origami pattern with concentrated hinges.The code
% will automatically generate compliant creases.

% thickness of two layers
tg=0.2*10^-6;
ts=0.80*10^-6;

% thickness of panel
tpanel=10*10^-6;

% width of creases
W=100*10^(-6);

% power input (mW)
qload=2.5;

% define the geometry
panelL=1*10^(-3);
undercut=W;

% generate a long stripe
ori.node0=[0 0 0; %1
      panelL 0 0;
      2*panelL+W 0 0;
      3*panelL+1.5*W 0 0;
      4*panelL+2*W 0 0;
      5*panelL+2.5*W 0 0;
      6*panelL+3*W 0 0; %7
      0 panelL 0; %8
      panelL panelL 0;
      2*panelL+W panelL 0;
      3*panelL+1.5*W panelL 0;
      4*panelL+2*W panelL 0;
      5*panelL+2.5*W panelL 0;
      6*panelL+3*W panelL 0; %14
      0*panelL,-0.6*panelL,-undercut;
      0*panelL,1.6*panelL,-undercut;
      8*panelL,-0.6*panelL,-undercut;
      8*panelL,1.6*panelL,-undercut;];
  
ori.panel0{1}=[1 2 9 8];
ori.panel0{2}=[2 3 10 9];
ori.panel0{3}=[3 4 11 10];
ori.panel0{4}=[4 5 12 11];
ori.panel0{5}=[5 6 13 12];
ori.panel0{6}=[6 7 14 13];

ori.panel0{7}=[15 16 18 17];


% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthMat=zeros(ori.oldCreaseNum,1);
ori.creaseWidthMat(3)=W;
ori.creaseWidthMat(6)=W;
ori.creaseWidthMat(9)=W;
ori.creaseWidthMat(12)=W;
ori.creaseWidthMat(15)=W;


% Compute the meshed geometry
ori.Mesh_CompliantCreaseGeometry()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=8*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;



%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickMat=[tpanel;tpanel;
                   tpanel;tpanel;
                   tpanel;tpanel;
                   tpanel;]; 
ori.panelW=W;


ori.creaseThickMat=zeros(ori.oldCreaseNum,1);
ori.creaseThickMat(3)=(tg+ts);
ori.creaseThickMat(6)=(tg+ts);
ori.creaseThickMat(9)=(tg+ts);
ori.creaseThickMat(12)=(tg+ts);
ori.creaseThickMat(15)=(tg+ts);


%% Assign Thermal Properties

ori.panelThermalConductMat = [0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;


% thickness of the surounding environment at RT
ori.t2RT=panelL/2; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerThermalLoading;

thermal.thermalStep=250;
thermal.tol=5*10^-7; 

thermal.supp=[1,1,1,1;
      2,1,1,1;
      3,1,1,1;
      4,1,1,1;
      25,1,1,1;
      26,1,1,1;
      27,1,1,1;
      28,1,1,1;];

thermal.thermalBoundaryPanelMat=[7];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
thermal.deltaAlpha(3)=50*10^(-6); 
thermal.deltaAlpha(6)=50*10^(-6); 
thermal.deltaAlpha(9)=50*10^(-6); 
thermal.deltaAlpha(12)=50*10^(-6); 
thermal.deltaAlpha(15)=50*10^(-6); 


thermal.Emat1=50*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=tg;
thermal.tmat2=ts;
thermal.videoOpen=1; % close the animation

% the target loading of crease heating
thermal.targetCreaseHeating=[3,qload/1000;
                             6,qload/1000;
                             9,qload/1000;
                             12,qload/1000;
                             15,qload/1000;];

ori.loadingController{1}={"ThermalLoading",thermal};


%% Solving the model
ori.Solver_Solve();

Uhis=thermal.Uhis;
U=Uhis(:,22,3)

