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


% thickness of two layer
t1=0.02*10^-3;
t2=0.02*10^-3;

% thickness of panel
tpanel=1*10^-3;

% width of crease
W=2*10^-3;

% undercut of XeF2 etching
underCut=2*10^-3;

% power input (W)
qload=0.06;

ori.node0=[-a -a 0;
           a -a 0;
           a a 0;
           -a a 0;
           -2*a -2*a 0;
           2*a -2*a 0;
           2*a 2*a 0;
           -2*a 2*a 0;
           -2*a 4*a 0;
           -a 5*a 0;
           a 5*a 0;
           2*a 4*a 0;
           -2*a 5*a 0;
           -a 6*a 0;
           a 6*a 0;
           2*a 5*a 0;
           -2*a -4*a 0;
           -a -5*a 0;
           a -5*a 0;
           2*a -4*a 0;
           -2*a -5*a 0;
           -a -6*a 0;
           a -6*a 0;
           2*a -5*a 0;
           -4*a -8*a -underCut;
           -4*a 8*a -underCut;
           5*a -8*a -underCut;
           5*a 8*a -underCut;
           ];
  
ori.panel0{1}=[1 2 3 4];
ori.panel0{2}=[5 1 4 8];
ori.panel0{3}=[2 6 7 3];
ori.panel0{4}=[3 7 12 11];
ori.panel0{5}=[4 3 11 10];
ori.panel0{6}=[8 4 10 9];
ori.panel0{7}=[9 10 14 13];
ori.panel0{8}=[10 11 15 14];
ori.panel0{9}=[11 12 16 15];
ori.panel0{10}=[5 17 18 1];
ori.panel0{11}=[18 19 2 1];
ori.panel0{12}=[19 20 6 2];
ori.panel0{13}=[21 22 18 17];
ori.panel0{14}=[22 23 19 18];
ori.panel0{15}=[23 24 20 19];
ori.panel0{16}=[25 27 28 26];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthMat=zeros(ori.oldCreaseNum,1);
ori.creaseWidthMat(1)=W;
ori.creaseWidthMat(2)=W;
ori.creaseWidthMat(3)=W;
ori.creaseWidthMat(4)=W;
ori.creaseWidthMat(6)=W;
ori.creaseWidthMat(7)=W;
ori.creaseWidthMat(8)=W;
ori.creaseWidthMat(10)=W;
ori.creaseWidthMat(11)=W;
ori.creaseWidthMat(14)=W;
ori.creaseWidthMat(13)=W;
ori.creaseWidthMat(15)=W;
ori.creaseWidthMat(17)=W;
ori.creaseWidthMat(19)=W;
ori.creaseWidthMat(21)=W;

ori.creaseWidthMat(27)=W;
ori.creaseWidthMat(29)=W;
ori.creaseWidthMat(26)=W;
ori.creaseWidthMat(28)=W;
ori.creaseWidthMat(30)=W;
ori.creaseWidthMat(34)=W;
ori.creaseWidthMat(36)=W;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=100*10^(-3); % plotting range
ori.displayRangeRatio=1; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

ori.panelThickMat=tpanel*ones(16,1); 
ori.panelW=W;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=1000;

ori.creaseThickMat=zeros(ori.oldCreaseNum,1);

ori.creaseThickMat(1)=t1+t2;
ori.creaseThickMat(2)=t1+t2;
ori.creaseThickMat(3)=t1+t2;
ori.creaseThickMat(4)=t1+t2;
ori.creaseThickMat(6)=t1+t2;
ori.creaseThickMat(7)=t1+t2;
ori.creaseThickMat(8)=t1+t2;
ori.creaseThickMat(10)=t1+t2;
ori.creaseThickMat(11)=t1+t2;
ori.creaseThickMat(14)=t1+t2;
ori.creaseThickMat(13)=t1+t2;
ori.creaseThickMat(15)=t1+t2;
ori.creaseThickMat(17)=t1+t2;
ori.creaseThickMat(19)=t1+t2;
ori.creaseThickMat(21)=t1+t2;

ori.creaseThickMat(27)=t1+t2;
ori.creaseThickMat(29)=t1+t2;
ori.creaseThickMat(26)=t1+t2;
ori.creaseThickMat(28)=t1+t2;
ori.creaseThickMat(30)=t1+t2;
ori.creaseThickMat(34)=t1+t2;
ori.creaseThickMat(36)=t1+t2;


%% Assign Thermal Properties

ori.panelThermalConductMat = 0.3*ones(16,1); 

ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=0.02; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerThermalLoading;

thermal.thermalStep=150;
thermal.tol=5*10^-5; 

thermal.supp=[1,1,1,1;
      2,1,1,1;
      3,1,1,1;
      4,1,1,1;
      61,1,1,1;
      62,1,1,1;
      63,1,1,1;
      64,1,1,1;];

thermal.thermalBoundaryPanelMat=[16];
thermal.roomTempNode=[];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
deltaAlpha=1000*10^(-6);
thermal.deltaAlpha(1)= deltaAlpha;
thermal.deltaAlpha(2)= deltaAlpha;
thermal.deltaAlpha(3)= deltaAlpha;
thermal.deltaAlpha(4)= deltaAlpha;
thermal.deltaAlpha(6)= deltaAlpha;
thermal.deltaAlpha(7)= deltaAlpha;
thermal.deltaAlpha(8)= deltaAlpha;
thermal.deltaAlpha(10)= deltaAlpha;

thermal.deltaAlpha(11)= deltaAlpha;
thermal.deltaAlpha(14)= deltaAlpha;
thermal.deltaAlpha(13)= deltaAlpha;
thermal.deltaAlpha(15)= deltaAlpha;
thermal.deltaAlpha(17)= deltaAlpha;
thermal.deltaAlpha(19)= deltaAlpha;
thermal.deltaAlpha(21)= deltaAlpha;


thermal.deltaAlpha(27)= deltaAlpha;
thermal.deltaAlpha(29)= deltaAlpha;
thermal.deltaAlpha(26)= deltaAlpha;
thermal.deltaAlpha(28)= deltaAlpha;
thermal.deltaAlpha(30)= deltaAlpha;
thermal.deltaAlpha(34)= deltaAlpha;
thermal.deltaAlpha(36)= deltaAlpha;


thermal.Emat1=2*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=t1;
thermal.tmat2=t2;
thermal.videoOpen=0;

% the target loading of crease heating
thermal.targetCreaseHeating=[
    1,qload;
    2,qload;
    3,qload;
    4,qload;
    6,qload;
    7,qload;
    8,qload;
    10,qload;
    
    11,qload;
    14,qload;
    13,qload;
    15,qload;
    17,qload;
    19,qload;
    21,qload;
    
    27,qload;
    29,qload;
    26,qload;
    28,qload;
    30,qload;
    34,qload;
    36,qload];

ori.loadingController{1}={"ThermalLoading",thermal};
ori.Solver_Solve();

Uhis_thermal1=thermal.Uhis;
x1=Uhis_thermal1(:,53,2) + ori.newNode(53,2);
x2=Uhis_thermal1(:,32,2) + ori.newNode(32,2);
distance_thermal1=(x1-x2).^2;


%% Insert another force loading
ori.continuingLoading=1;

nr=ControllerNRLoading;

nr.increStep=50;
nr.tol=10^-6;
nr.iterMax=50;

nr.supp=[1,1,1,1;
      2,1,1,1;
      3,1,1,1;
      4,1,1,1;
      61,1,1,1;
      62,1,1,1;
      63,1,1,1;
      64,1,1,1;];

loadMag=0.9*(a*a*tpanel*2000*9.8)/5;
nr.load=[53,0,-loadMag,0;
      54,0,-loadMag,0;
      31,0,loadMag,0;
      32,0,loadMag,0;];
  
nr.videoOpen=0;

ori.loadingController{1}={"NR",nr};
ori.Solver_Solve();


Uhis_NR=nr.Uhis;
x1=Uhis_NR(:,53,2) + ori.newNode(53,2);
x2=Uhis_NR(:,32,2) + ori.newNode(32,2);
distance_NR=(x1-x2).^2;

%% The second thermal loading step
ori.continuingLoading=1;
thermal.thermalStep=120;

qload=qload;
thermal.targetCreaseHeating=[
    1,qload;
    2,qload;
    3,qload;
    4,qload;
    6,qload;
    7,qload;
    8,qload;
    10,qload;
    
    11,qload;
    14,qload;
    13,qload;
    15,qload;
    17,qload;
    19,qload;
    21,qload;
    
    27,qload;
    29,qload;
    26,qload;
    28,qload;
    30,qload;
    34,qload;
    36,qload];

ori.loadingController{1}={"ThermalLoading",thermal};
ori.Solver_Solve();

Uhis_thermal2=thermal.Uhis;
x1=Uhis_thermal2(:,53,2) + ori.newNode(53,2);
x2=Uhis_thermal2(:,32,2) + ori.newNode(32,2);
distance_thermal2=(x1-x2).^2;

toc
