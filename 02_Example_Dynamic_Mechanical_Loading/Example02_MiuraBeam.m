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
clear all;clc;close all;
ori=OrigamiSolver;
tic

%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

a=50*10^(-3);
b=50*10^(-3);
gama=80*pi/180;
m=4;
n=2;
Ext=0.6;
[node0,panel0]=GenerateMiuraSheet(a,b,gama,m,n,Ext);

ori.node0=node0;
ori.panel0=panel0;

% Ploting option
ori.viewAngle1=15;
ori.viewAngle2=30;
ori.displayRange=220*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();
ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;

%% Meshing of the origami model
% Define the crease width
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=6*10^(-3);
ori.creaseWidthVec(4)=6*10^(-3);
ori.creaseWidthVec(6)=6*10^(-3);

ori.creaseWidthVec(9)=6*10^(-3);
ori.creaseWidthVec(10)=6*10^(-3);
ori.creaseWidthVec(11)=6*10^(-3);

ori.creaseWidthVec(15)=6*10^(-3);
ori.creaseWidthVec(14)=6*10^(-3);
ori.creaseWidthVec(16)=6*10^(-3);
ori.creaseWidthVec(20)=6*10^(-3);

% Define other parameters for meshing
ori.mesh2D3D=3; % load a 3D meshing method
ori.compliantCreaseOpen=0;

% Compute the meshed geometry
ori.Mesh_Mesh()
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;

%% Assign Mechanical Properties

ori.panelE=1000*10^6; 
ori.creaseE=2000*10^6; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3;
ori.panelThickVec=[1;1;1;1;...
                   1;1;1;1;]*1000*10^(-6); 
               
ori.panelW=6*10^-3;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=100*10^(-6);
ori.creaseThickVec(4)=100*10^(-6);
ori.creaseThickVec(6)=100*10^(-6);

ori.creaseThickVec(9)=100*10^(-6);
ori.creaseThickVec(10)=100*10^(-6);
ori.creaseThickVec(11)=100*10^(-6);

ori.creaseThickVec(15)=100*10^(-6);
ori.creaseThickVec(14)=100*10^(-6);
ori.creaseThickVec(16)=100*10^(-6);
ori.creaseThickVec(20)=100*10^(-6);

%% Setup the loading controller
rhoSU8=1200;
ori.densityCrease=rhoSU8;
ori.densityPanel=rhoSU8;

dynamics=ControllerDynamics();

dynamics.supp=[1,1,1,1;
      3,0,0,1;
      4,0,0,1;]; 
dynamics.dt=10^-5;
  
step=40000;
TimeVec=(1:step)*10^-5;
dynamics.Fext=zeros(step,23,3);
dynamics.rotTargetAngle=pi*ones(step,22);
% This will force the system to spring back to flat

% Apply a step force
% dynamics.Fext(:,6,3)=0.0000001;
% dynamics.Fext(:,7,3)=0.0000001;
% Apply a step change in stress free folding angle
% dynamics.rotTargetAngle(:,3)=pi+pi/4;

% set up the damping coefficient
dynamics.alpha=0.003;
dynamics.beta=0.003;

% ploting option
dynamics.plotOpen=0;
dynamics.videoOpen=1;
dynamics.videoCropRate=100;

% Solve the solution
ori.loadingController{1}={"Dynamics",dynamics};
ori.Solver_Solve()
toc
