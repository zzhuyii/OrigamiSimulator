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

a=20*10^(-3);
b=4*10^(-3);
c=10*10^(-3);
z=3*10^(-3);

ori.node0=[b 0 0;
      a+b 0 0;
      2*a+b 0 0;
      3*a+b 0 0;
      4*a+b 0 0; 
      5*a+b 0 0; %6
      
      0 c 0;
      a c 0;
      2*a c 0;
      3*a+2*b c 0;
      4*a+2*b c 0; 
      5*a+2*b c 0; %12
      
      b 2*c 0;
      a+b 2*c 0;
      2*a+b 2*c 0;
      3*a+b 2*c 0;
      4*a+b 2*c 0; 
      5*a+b 2*c 0; %18
      
      b 0 z;
      -z+1*a+b 0 z;
      -z+2*a+b 0 z;
      z+3*a+b 0 z;
      z+4*a+b 0 z; 
      z+5*a+b 0 z; %6
      
      0 c z;
      -z+1*a c z;
      -z+2*a c z;
      z+3*a+2*b c z;
      z+4*a+2*b c z; 
      z+5*a+2*b c z; %12
      
      b 2*c z;
      -z+1*a+b 2*c z;
      -z+2*a+b 2*c z;
      z+3*a+b 2*c z;
      z+4*a+b 2*c z;
      z+5*a+b 2*c z;];
  
ori.panel0{1}=[1 2 8 7];
ori.panel0{2}=[2 3 9 8];
ori.panel0{3}=[3 4 10 9];
ori.panel0{4}=[4 5 11 10];
ori.panel0{5}=[5 6 12 11];
ori.panel0{6}=[7 8 14 13];
ori.panel0{7}=[8 9 15 14];
ori.panel0{8}=[9 10 16 15];
ori.panel0{9}=[10 11 17 16];
ori.panel0{10}=[11 12 18 17];

ori.panel0{11}=[1 2 8 7]+[18 18 18 18];
ori.panel0{12}=[2 3 9 8]+[18 18 18 18];
ori.panel0{13}=[3 4 10 9]+[18 18 18 18];
ori.panel0{14}=[4 5 11 10]+[18 18 18 18];
ori.panel0{15}=[5 6 12 11]+[18 18 18 18];
ori.panel0{16}=[7 8 14 13]+[18 18 18 18];
ori.panel0{17}=[8 9 15 14]+[18 18 18 18];
ori.panel0{18}=[9 10 16 15]+[18 18 18 18];
ori.panel0{19}=[10 11 17 16]+[18 18 18 18];
ori.panel0{20}=[11 12 18 17]+[18 18 18 18];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);

ori.creaseWidthVec(4)=3*10^(-3);
ori.creaseWidthVec(7)=3*10^(-3);
ori.creaseWidthVec(10)=3*10^(-3);
ori.creaseWidthVec(13)=3*10^(-3);
ori.creaseWidthVec(16)=3*10^(-3);

ori.creaseWidthVec(3)=3*10^(-3);
ori.creaseWidthVec(18)=3*10^(-3);
ori.creaseWidthVec(6)=3*10^(-3);
ori.creaseWidthVec(20)=3*10^(-3);
ori.creaseWidthVec(9)=3*10^(-3);
ori.creaseWidthVec(22)=3*10^(-3);
ori.creaseWidthVec(12)=3*10^(-3);
ori.creaseWidthVec(24)=3*10^(-3);

ori.creaseWidthVec(4+27)=3*10^(-3);
ori.creaseWidthVec(7+27)=3*10^(-3);
ori.creaseWidthVec(10+27)=3*10^(-3);
ori.creaseWidthVec(13+27)=3*10^(-3);
ori.creaseWidthVec(16+27)=3*10^(-3);

ori.creaseWidthVec(3+27)=3*10^(-3);
ori.creaseWidthVec(18+27)=3*10^(-3);
ori.creaseWidthVec(6+27)=3*10^(-3);
ori.creaseWidthVec(20+27)=3*10^(-3);
ori.creaseWidthVec(9+27)=3*10^(-3);
ori.creaseWidthVec(22+27)=3*10^(-3);
ori.creaseWidthVec(12+27)=3*10^(-3);
ori.creaseWidthVec(24+27)=3*10^(-3);

% close the compliant crease
ori.compliantCreaseOpen=0;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=15;
ori.viewAngle2=15;
ori.displayRange=70*10^(-3); % plotting range
ori.displayRangeRatio=0.1; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=ones(20,1)*500*10^(-6);
ori.panelW=3*10^-3;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);

ori.creaseThickVec(4)=100*10^(-6);
ori.creaseThickVec(7)=100*10^(-6);
ori.creaseThickVec(10)=100*10^(-6);
ori.creaseThickVec(13)=100*10^(-6);
ori.creaseThickVec(16)=100*10^(-6);

ori.creaseThickVec(3)=100*10^(-6);
ori.creaseThickVec(18)=100*10^(-6);
ori.creaseThickVec(6)=100*10^(-6);
ori.creaseThickVec(20)=100*10^(-6);
ori.creaseThickVec(9)=100*10^(-6);
ori.creaseThickVec(22)=100*10^(-6);
ori.creaseThickVec(12)=100*10^(-6);
ori.creaseThickVec(24)=100*10^(-6);

ori.creaseThickVec(4+27)=100*10^(-6);
ori.creaseThickVec(7+27)=100*10^(-6);
ori.creaseThickVec(10+27)=100*10^(-6);
ori.creaseThickVec(13+27)=100*10^(-6);
ori.creaseThickVec(16+27)=100*10^(-6);

ori.creaseThickVec(3+27)=100*10^(-6);
ori.creaseThickVec(18+27)=100*10^(-6);
ori.creaseThickVec(6+27)=100*10^(-6);
ori.creaseThickVec(20+27)=100*10^(-6);
ori.creaseThickVec(9+27)=100*10^(-6);
ori.creaseThickVec(22+27)=100*10^(-6);
ori.creaseThickVec(12+27)=100*10^(-6);
ori.creaseThickVec(24+27)=100*10^(-6);


%% setup panel contact information
ori.contactOpen=1;
ori.ke=0.002;
ori.d0edge=1*10^-3;
ori.d0center=1*10^-3;


%% setup the loading process

% define the self folding step
selfFold=ControllerSelfFolding;

% Assign zero strain position for creases during self-folding
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);

ratio=0.5;
selfFold.targetRotZeroStrain(4)=pi+ratio*pi;
selfFold.targetRotZeroStrain(7)=pi-ratio*pi;
selfFold.targetRotZeroStrain(10)=pi+ratio*pi;
selfFold.targetRotZeroStrain(13)=pi-ratio*pi;
selfFold.targetRotZeroStrain(16)=pi+ratio*pi;

selfFold.targetRotZeroStrain(3)=pi+ratio*pi;
selfFold.targetRotZeroStrain(18)=pi+ratio*pi;
selfFold.targetRotZeroStrain(6)=pi-ratio*pi;
selfFold.targetRotZeroStrain(20)=pi-ratio*pi;
selfFold.targetRotZeroStrain(9)=pi-ratio*pi;
selfFold.targetRotZeroStrain(22)=pi-ratio*pi;
selfFold.targetRotZeroStrain(12)=pi+ratio*pi;
selfFold.targetRotZeroStrain(24)=pi+ratio*pi;

selfFold.targetRotZeroStrain(4+27)=pi+ratio*pi;
selfFold.targetRotZeroStrain(7+27)=pi-ratio*pi;
selfFold.targetRotZeroStrain(10+27)=pi+ratio*pi;
selfFold.targetRotZeroStrain(13+27)=pi-ratio*pi;
selfFold.targetRotZeroStrain(16+27)=pi+ratio*pi;

selfFold.targetRotZeroStrain(3+27)=pi+ratio*pi;
selfFold.targetRotZeroStrain(18+27)=pi+ratio*pi;
selfFold.targetRotZeroStrain(6+27)=pi-ratio*pi;
selfFold.targetRotZeroStrain(20+27)=pi-ratio*pi;
selfFold.targetRotZeroStrain(9+27)=pi-ratio*pi;
selfFold.targetRotZeroStrain(22+27)=pi-ratio*pi;
selfFold.targetRotZeroStrain(12+27)=pi+ratio*pi;
selfFold.targetRotZeroStrain(24+27)=pi+ratio*pi;

selfFold.supp=[1,1,1,1;
          2,1,1,1;
          13,0,0,1;
          14,0,0,1;
          19,1,1,1;
          20,1,1,1;
          31,0,0,1;
          32,0,0,1;];

selfFold.increStep=30;
selfFold.tol=10^-6;
selfFold.iterMax=50;      



% define a NR loading step
nr=ControllerNRLoading;

nr.supp=[2,1,1,1;
      14,1,1,1;
      5,1,1,1;
      17,1,1,1;
      32,1,1,0;
      20,1,1,0;
      23,1,1,0;
      35,1,1,0;];
      
nr.nonRigidSupport=1;
ksup=5;
nr.suppElastic=[19,3,ksup;
             20,3,ksup;
             31,3,ksup;
             32,3,ksup;
             23,3,ksup;
             24,3,ksup;
             35,3,ksup;
             36,3,ksup;];
      
loadForce=5*10^(-3);
nr.load=[21,0,0,loadForce;
      22,0,0,loadForce;
      33,0,0,loadForce;
      34,0,0,loadForce;]; 
     
nr.increStep=30;
nr.tol=10^-6;
nr.iterMax=50;


ori.loadingController{1}={"SelfFold",selfFold};
ori.loadingController{2}={"NR",nr};

%% Solving the model
ori.Solver_Solve();


