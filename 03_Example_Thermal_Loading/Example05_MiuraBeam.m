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
W=200*10^(-6);

% power input (mW)
qload=9;

% generate a miura beam
a=1*10^(-3);
b=1*10^(-3);
gama=80*pi/180;
m=2;
n=5;
Ext=1;

[ori.node0,ori.panel0]=GenerateMiuraSheet(a,b,gama,m,n,Ext);
ori.node0(19,:)=[-8,-8,-0.1]*10^(-3);
ori.node0(20,:)=[-8,8,-0.1]*10^(-3);
ori.node0(21,:)=[8,8,-0.1]*10^(-3);
ori.node0(22,:)=[8,-8,-0.1]*10^(-3);

ori.panel0{11}=[19,20,21,22];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;
ori.creaseWidthVec(4)=W;
ori.creaseWidthVec(19)=W;
ori.creaseWidthVec(6)=W;
ori.creaseWidthVec(7)=W;
ori.creaseWidthVec(21)=W;
ori.creaseWidthVec(9)=W;
ori.creaseWidthVec(10)=W;
ori.creaseWidthVec(23)=W;
ori.creaseWidthVec(12)=W;
ori.creaseWidthVec(13)=W;
ori.creaseWidthVec(25)=W;
ori.creaseWidthVec(15)=W;


% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=6*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;



%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[tpanel;tpanel;
                   tpanel;tpanel;
                   tpanel;tpanel;
                   tpanel;tpanel;
                   tpanel;tpanel;tpanel;]; 
ori.panelW=W;


ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(tg+ts);
ori.creaseThickVec(4)=(tg+ts);
ori.creaseThickVec(19)=(tg+ts);
ori.creaseThickVec(6)=(tg+ts);
ori.creaseThickVec(7)=(tg+ts);
ori.creaseThickVec(21)=(tg+ts);
ori.creaseThickVec(9)=(tg+ts);
ori.creaseThickVec(10)=(tg+ts);
ori.creaseThickVec(23)=(tg+ts);
ori.creaseThickVec(12)=(tg+ts);
ori.creaseThickVec(13)=(tg+ts);
ori.creaseThickVec(25)=(tg+ts);
ori.creaseThickVec(15)=(tg+ts);


%% Assign Thermal Properties

ori.panelThermalConductVec = [0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;0.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;


% thickness of the surounding environment at RT
ori.t2RT=(a+b)/2; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=250;
thermal.tol=5*10^-7; 

thermal.supp=[1,1,1,1;
      4,1,0,1;
      22,0,0,1;
      23,0,0,1;
      41,1,1,1;
      42,1,1,1;
      43,1,1,1;
      44,1,1,1];

thermal.thermalBoundaryPanelVec=[11];
thermal.roomTempNode=[];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
thermal.deltaAlpha(3)=-50*10^(-6); 
thermal.deltaAlpha(6)=50*10^(-6); 
thermal.deltaAlpha(4)=50*10^(-6); 
thermal.deltaAlpha(19)=50*10^(-6); 

thermal.deltaAlpha(7)=-50*10^(-6); 
thermal.deltaAlpha(21)=-50*10^(-6); 
thermal.deltaAlpha(9)=-50*10^(-6); 

thermal.deltaAlpha(10)=50*10^(-6); 
thermal.deltaAlpha(23)=50*10^(-6); 
thermal.deltaAlpha(12)=50*10^(-6); 

thermal.deltaAlpha(13)=-50*10^(-6); 
thermal.deltaAlpha(25)=-50*10^(-6); 
thermal.deltaAlpha(15)=-50*10^(-6); 


thermal.Emat1=50*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=tg;
thermal.tmat2=ts;
thermal.videoOpen=0; % close the animation

% the target loading of crease heating
thermal.targetCreaseHeating=[3,qload/1000;
                             4,qload/1000;
                             6,qload/1000;
                             19,qload/1000;
                             7,qload/1000;
                             21,qload/1000;
                             9,qload/1000;
                             10,qload/1000;
                             23,qload/1000;
                             12,qload/1000;
                             13,qload/1000;
                             25,qload/1000;
                             15,qload/1000;];

ori.loadingController{1}={"ElectroThermal",thermal};


%% Solving the model
ori.Solver_Solve();

Uhis=thermal.Uhis;
U=Uhis(:,19,3)
toc