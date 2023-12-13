%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer: Yi Zhu
% Advisor: Evgueni T. Filipov
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% origami simulator presented in this package. 
%
% Reference:
% [1] Yi Zhu, Evgueni T. Filipov (2021). 'Sequentially Working Origami 
%     Multi-Physics Simulator (SWOMPS): A Versatile Implementation',
%     ASME IDETC-CIE Conference. DETC2021-68042. 
% [2] Y. Zhu, E. T. Filipov (2021). 'Rapid Multi-Physic Simulation for 
%     Electro-Thermal Origami Robotic Systems'  International Journal of 
%     Mechanical Sciences, 202-203, 106537.
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
%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%

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
W=200*10^(-6);

% power input (mW)
qload=9;

% generate a miura beam
a=1*10^(-3);
b=1*10^(-3);
gama=60*pi/180;
m=2;
n=2;
Ext=1;

[ori.node0,ori.panel0]=GenerateMiuraSheet(a,b,gama,m,n,Ext);

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;
ori.creaseWidthVec(4)=W;
ori.creaseWidthVec(10)=W;
ori.creaseWidthVec(6)=W;


% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=280;
ori.viewAngle2=30;
ori.displayRange=3*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;



%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[tpanel;tpanel;tpanel;tpanel;]; 
ori.panelW=W;


ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(tg+ts);
ori.creaseThickVec(6)=(tg+ts);
ori.creaseThickVec(4)=(tg+ts);
ori.creaseThickVec(10)=(tg+ts);


%% Assign Thermal Properties

ori.panelThermalConductVec = [0.3;0.3;0.3;0.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;


% thickness of the surounding environment at RT
ori.t2RT=(a+b)/2; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=250;
thermal.tol=5*10^-7; 

thermal.supp=[5,1,1,1;
      8,1,0,1;
      14,0,0,1;
      15,0,0,1;];


thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
dAlpha=50*10^(-6);
thermal.deltaAlpha(3)=-dAlpha; 
thermal.deltaAlpha(6)=dAlpha; 
thermal.deltaAlpha(4)=dAlpha; 
thermal.deltaAlpha(10)=dAlpha; 

thermal.Emat1=50*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=tg;
thermal.tmat2=ts;
thermal.videoOpen=0; % close the animation

% the target loading of crease heating
thermal.targetCreaseHeating=[3,qload/1000;
                             4,qload/1000;
                             6,qload/1000;
                             10,qload/1000];

ori.loadingController{1}={"ElectroThermal",thermal};


%% Solving the model
ori.Solver_Solve();
