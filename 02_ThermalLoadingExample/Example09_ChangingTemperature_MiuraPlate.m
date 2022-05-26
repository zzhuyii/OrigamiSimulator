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
tpanel=5*10^-6;

% width of creases
W=200*10^(-6);

% generate a miura beam
a=1*10^(-3);
b=1*10^(-3);
gama=60*pi/180;
m=6;
n=6;
Ext=1.0;

[ori.node0,ori.panel0]=GenerateMiuraSheet(a,b,gama,m,n,Ext);

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model


% Find the crease numbers for those internal folding crease
A=ori.oldCreaseType;
internalCreaseNum=find(A==2);

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(internalCreaseNum)=W;


% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=0;
ori.viewAngle2=90;

ori.viewAngle1=60;
ori.viewAngle2=40;
ori.width=1200;
ori.height=1200;



ori.displayRange=8*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

% ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
% ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;



%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

% find the number of panels in the original origami pattern
A=size(ori.panel0);
panelNum=A(2);

% set up the thickness of origami panels
ori.panelThickVec=tpanel*ones(panelNum,1);
ori.panelW=W;

% Define the thickness of the creases
ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(internalCreaseNum)=(tg+ts);



%% Assign Thermal Properties

ori.panelThermalConductVec = 0.3*ones(panelNum,1); 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;


% thickness of the surounding environment at RT
ori.t2RT=(a+b)/2; 



%% Setup the loading controller
% Run one initialization step to help convergence
temperatureIni=ControllerChangingTemperature;

temperatureIni.thermalStep=2;
temperatureIni.tol=5*10^-7; 

temperatureIni.supp=[340,1,1,1;
      339,1,0,1;
      346,0,0,1;
      345,0,0,1;];


temperatureIni.deltaAlpha=zeros(ori.oldCreaseNum,1);
dAlpha=20*10^(-6);

temperatureIni.deltaAlpha(internalCreaseNum)=dAlpha; 

temperatureIni.deltaAlpha(15)=-dAlpha; 
temperatureIni.deltaAlpha(31)=-dAlpha; 
temperatureIni.deltaAlpha(42)=-dAlpha; 
temperatureIni.deltaAlpha(57)=-dAlpha; 
temperatureIni.deltaAlpha(68)=-dAlpha; 

temperatureIni.deltaAlpha(9)=-dAlpha; 
temperatureIni.deltaAlpha(27)=-dAlpha; 
temperatureIni.deltaAlpha(38)=-dAlpha; 
temperatureIni.deltaAlpha(53)=-dAlpha; 
temperatureIni.deltaAlpha(64)=-dAlpha; 

temperatureIni.deltaAlpha(3)=-dAlpha; 
temperatureIni.deltaAlpha(23)=-dAlpha; 
temperatureIni.deltaAlpha(34)=-dAlpha; 
temperatureIni.deltaAlpha(49)=-dAlpha; 
temperatureIni.deltaAlpha(60)=-dAlpha; 



temperatureIni.deltaAlpha(13)=-dAlpha; 
temperatureIni.deltaAlpha(28)=-dAlpha; 
temperatureIni.deltaAlpha(41)=-dAlpha; 
temperatureIni.deltaAlpha(54)=-dAlpha; 
temperatureIni.deltaAlpha(67)=-dAlpha; 
temperatureIni.deltaAlpha(80)=-dAlpha; 

temperatureIni.deltaAlpha(7)=-dAlpha; 
temperatureIni.deltaAlpha(24)=-dAlpha; 
temperatureIni.deltaAlpha(37)=-dAlpha; 
temperatureIni.deltaAlpha(50)=-dAlpha; 
temperatureIni.deltaAlpha(63)=-dAlpha; 
temperatureIni.deltaAlpha(76)=-dAlpha; 


temperatureIni.Emat1=50*10^9; 
temperatureIni.Emat2=2*10^9;
temperatureIni.tmat1=tg;
temperatureIni.tmat2=ts;
temperatureIni.plotOpen=0;
temperatureIni.videoOpen=0; % close the animation

% the target loading of crease heating
temperatureIni.targetAmbientTemperature=ori.RT;


% Change the ambient temperature to 150C
temperature=ControllerChangingTemperature;

temperature.thermalStep=100; % 1 step increase 3C temperature
temperature.tol=5*10^-7; 

temperature.supp=[340,1,1,1;
      339,1,0,1;
      346,0,0,1;
      345,0,0,1;];


temperature.deltaAlpha=zeros(ori.oldCreaseNum,1);
temperature.deltaAlpha(internalCreaseNum)=dAlpha; 

temperature.deltaAlpha(15)=-dAlpha; 
temperature.deltaAlpha(31)=-dAlpha; 
temperature.deltaAlpha(42)=-dAlpha; 
temperature.deltaAlpha(57)=-dAlpha; 
temperature.deltaAlpha(68)=-dAlpha; 

temperature.deltaAlpha(9)=-dAlpha; 
temperature.deltaAlpha(27)=-dAlpha; 
temperature.deltaAlpha(38)=-dAlpha; 
temperature.deltaAlpha(53)=-dAlpha; 
temperature.deltaAlpha(64)=-dAlpha; 

temperature.deltaAlpha(3)=-dAlpha; 
temperature.deltaAlpha(23)=-dAlpha; 
temperature.deltaAlpha(34)=-dAlpha; 
temperature.deltaAlpha(49)=-dAlpha; 
temperature.deltaAlpha(60)=-dAlpha; 



temperature.deltaAlpha(13)=-dAlpha; 
temperature.deltaAlpha(28)=-dAlpha; 
temperature.deltaAlpha(41)=-dAlpha; 
temperature.deltaAlpha(54)=-dAlpha; 
temperature.deltaAlpha(67)=-dAlpha; 
temperature.deltaAlpha(80)=-dAlpha; 

temperature.deltaAlpha(7)=-dAlpha; 
temperature.deltaAlpha(24)=-dAlpha; 
temperature.deltaAlpha(37)=-dAlpha; 
temperature.deltaAlpha(50)=-dAlpha; 
temperature.deltaAlpha(63)=-dAlpha; 
temperature.deltaAlpha(76)=-dAlpha; 


temperature.Emat1=50*10^9; 
temperature.Emat2=2*10^9;
temperature.tmat1=tg;
temperature.tmat2=ts;
temperature.plotOpen=0;
temperature.videoOpen=0; % close the animation

% the target loading of crease heating
temperature.targetAmbientTemperature=3*temperature.thermalStep;

                         
% set up the four step loading           
ori.loadingController{1}={"ChangingTemperature",temperatureIni};
ori.loadingController{2}={"ChangingTemperature",temperature};

%% Solving the model
ori.Solver_Solve();
ori.Plot_DeformedShape(ori.newNode,ori.newNode+ori.currentU)
