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

tic
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
qload=4;

% generate crane pattern
a=1*10^(-3);
Node=[-3.2*a 0 0;
      -0.7*a 0 0;
      0 0 0;
      0.8*a 0 0;
      2.3*a 0 0;
      3.4*a 0 0; % 6 end of center
      -1.1*a a 0;
      1.1*a a 0;
      2.6*a 0.5*a 0;
      0 4*a 0; % 10 end of top
      -1.1*a -a 0;
      1.1*a -a 0;
      2.6*a -0.5*a 0;
      0*a -4*a 0; % 10 end of top
      ]; 
      
Panel{1}=[1 11 2];
Panel{2}=[1 2 7];
Panel{3}=[2 3 7];
Panel{4}=[2 11 3];
Panel{5}=[11 12 3];
Panel{6}=[7 3 8];
Panel{7}=[3 4 8];
Panel{8}=[3 12 4];
Panel{9}=[12 13 5 4];
Panel{10}=[4 5 9 8];
Panel{11}=[5 13 6];
Panel{12}=[5 6 9];
Panel{13}=[7 8 10];
Panel{14}=[11 14 12];

ori.node0=Node;
ori.panel0=Panel;

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(1)=W;
ori.creaseWidthVec(3)=W;
ori.creaseWidthVec(5)=W;
ori.creaseWidthVec(6)=W;
ori.creaseWidthVec(7)=W;
ori.creaseWidthVec(8)=W;
ori.creaseWidthVec(9)=W;
ori.creaseWidthVec(11)=W;
ori.creaseWidthVec(10)=W;
ori.creaseWidthVec(12)=W;
ori.creaseWidthVec(13)=W;
ori.creaseWidthVec(14)=W;
ori.creaseWidthVec(15)=W;
ori.creaseWidthVec(18)=W;
ori.creaseWidthVec(17)=W;
ori.creaseWidthVec(19)=W;
ori.creaseWidthVec(21)=W;


% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=10;
ori.viewAngle2=30;
ori.displayRange=4*10^(-3); % plotting range
ori.displayRangeRatio=1; % plotting range in the negative axis

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
                   tpanel;tpanel;
                   tpanel;tpanel;
                   tpanel;tpanel]; 
ori.panelW=W;


ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(1)=(tg+ts);
ori.creaseThickVec(3)=(tg+ts);
ori.creaseThickVec(5)=(tg+ts);
ori.creaseThickVec(6)=(tg+ts);
ori.creaseThickVec(7)=(tg+ts);
ori.creaseThickVec(8)=(tg+ts);
ori.creaseThickVec(9)=(tg+ts);
ori.creaseThickVec(11)=(tg+ts);
ori.creaseThickVec(10)=(tg+ts);
ori.creaseThickVec(12)=(tg+ts);
ori.creaseThickVec(13)=(tg+ts);
ori.creaseThickVec(14)=(tg+ts);
ori.creaseThickVec(15)=(tg+ts);
ori.creaseThickVec(18)=(tg+ts);
ori.creaseThickVec(17)=(tg+ts);
ori.creaseThickVec(19)=(tg+ts);
ori.creaseThickVec(21)=(tg+ts);


%% Assign Thermal Properties

ori.panelThermalConductVec = [0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;
                              0.3;0.3;]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;


% thickness of the surounding environment at RT
ori.t2RT=a; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=500;
thermal.tol=5*10^-7; 

thermal.supp=[13,1,1,1;
      14,0,1,1;
      16,0,0,1;
      18,0,0,1;];

thermal.thermalBoundaryPanelVec=[];
thermal.roomTempNode=[];

dAlpha=200*10^(-6);
thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
thermal.deltaAlpha(1)=-dAlpha; 
thermal.deltaAlpha(3)=dAlpha; 
thermal.deltaAlpha(5)=dAlpha; 
thermal.deltaAlpha(6)=dAlpha; 

thermal.deltaAlpha(7)=-dAlpha;
thermal.deltaAlpha(8)=-dAlpha;
thermal.deltaAlpha(10)=-dAlpha;
thermal.deltaAlpha(12)=-dAlpha;

thermal.deltaAlpha(13)=dAlpha; 
thermal.deltaAlpha(14)=dAlpha; 
thermal.deltaAlpha(15)=dAlpha; 
thermal.deltaAlpha(18)=-dAlpha; 

thermal.deltaAlpha(17)=-dAlpha; 
thermal.deltaAlpha(19)=-dAlpha;
thermal.deltaAlpha(21)=dAlpha; 

thermal.deltaAlpha(9)=dAlpha;
thermal.deltaAlpha(11)=dAlpha;


thermal.Emat1=50*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=tg;
thermal.tmat2=ts;
thermal.videoOpen=0; % close the animation

% the target loading of crease heating
thermal.targetCreaseHeating=[
                             6,qload/1000;
                             13,qload/1000;
                             
                             1,4*qload/1000;
                             18,4*qload/1000;
                             
                             3,qload/1000*1.5;
                             5,qload/1000*1.5;
                             14,qload/1000*1.5;
                             15,qload/1000*1.5;
                             
                             7,qload/1000*2;
                             8,qload/1000*2;
                             10,qload/1000*2;
                             12,qload/1000*2;
                                                          
                             17,qload/1000/2;
                             19,qload/1000/2;
                             
                             21,qload/1000;
                             
                             9,qload/1000;
                             11,qload/1000;];
                         

% When deactivate certain creses
% thermal.targetCreaseHeating=[
%                              6,qload/1000;
%                              13,qload/1000;
%                              
%                              1,4*qload/1000;
%                              18,4*qload/1000;
%                              
%                              3,qload/1000*1.5;
%                              5,qload/1000*1.5;
%                              14,qload/1000*1.5;
%                              15,qload/1000*1.5;
%                              
%                              7,qload/1000*2;
%                              %8,qload/1000*2;
%                              10,qload/1000*2;
%                              12,qload/1000*2;
%                                                           
%                              %17,qload/1000/2;
%                              19,qload/1000/2;
%                              
%                              21,qload/1000;
%                              
%                              9,qload/1000/2;
%                              11,qload/1000/2;];

ori.loadingController{1}={"ElectroThermal",thermal};


%% Solving the model
ori.Solver_Solve();
toc
