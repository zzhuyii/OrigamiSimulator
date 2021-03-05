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

tic
a=12*10^(-3);
b=30*10^(-3);
c=45*10^(-3);
ori.node0=[0 0 0;
      a*cos(0) a*sin(0) 0;
      a*cos(0.25*pi) a*sin(0.25*pi) 0;
      a*cos(0.5*pi) a*sin(0.5*pi) 0;
      a*cos(0.75*pi) a*sin(0.75*pi) 0;
      a*cos(pi) a*sin(pi) 0;
      a*cos(1.25*pi) a*sin(1.25*pi) 0;
      a*cos(1.5*pi) a*sin(1.5*pi) 0;
      a*cos(1.75*pi) a*sin(1.75*pi) 0;% end of inner loop
      b*cos((1/16)*pi) b*sin((1/16)*pi) 0;
      b*cos((0.125+1/16)*pi) b*sin((0.125+1/16)*pi) 0;
      b*cos((0.25+1/16)*pi) b*sin((0.25+1/16)*pi) 0;
      b*cos((0.375+1/16)*pi) b*sin((0.375+1/16)*pi) 0;
      b*cos((0.5+1/16)*pi) b*sin((0.5+1/16)*pi) 0;
      b*cos((0.625+1/16)*pi) b*sin((0.625+1/16)*pi) 0;
      b*cos((0.75+1/16)*pi) b*sin((0.75+1/16)*pi) 0;
      b*cos((0.875+1/16)*pi) b*sin((0.875+1/16)*pi) 0;
      b*cos((1+1/16)*pi) b*sin((1+1/16)*pi) 0;
      b*cos((1.125+1/16)*pi) b*sin((1.125+1/16)*pi) 0;
      b*cos((1.25+1/16)*pi) b*sin((1.25+1/16)*pi) 0;
      b*cos((1.375+1/16)*pi) b*sin((1.375+1/16)*pi) 0;
      b*cos((1.5+1/16)*pi) b*sin((1.5+1/16)*pi) 0;
      b*cos((1.625+1/16)*pi) b*sin((1.625+1/16)*pi) 0;
      b*cos((1.75+1/16)*pi) b*sin((1.75+1/16)*pi) 0;
      b*cos((1.875+1/16)*pi) b*sin((1.875+1/16)*pi) 0; % end of first outer loop
      b*cos(0) b*sin(0) 0;
      b*cos(0.25*pi) b*sin(0.25*pi) 0;
      b*cos(0.5*pi) b*sin(0.5*pi) 0;
      b*cos(0.75*pi) b*sin(0.75*pi) 0;
      b*cos(pi) b*sin(pi) 0;
      b*cos(1.25*pi) b*sin(1.25*pi) 0;
      b*cos(1.5*pi) b*sin(1.5*pi) 0;
      b*cos(1.75*pi) b*sin(1.75*pi) 0; % end of second outer loop
      c*cos(1/8*pi) c*sin(1/8*pi) 0;
      c*cos((0.25+1/8)*pi) c*sin((0.25+1/8)*pi) 0;
      c*cos((0.5+1/8)*pi) c*sin((0.5+1/8)*pi) 0;
      c*cos((0.75+1/8)*pi) c*sin((0.75+1/8)*pi) 0;
      c*cos((1+1/8)*pi) c*sin((1+1/8)*pi) 0;
      c*cos((1.25+1/8)*pi) c*sin((1.25+1/8)*pi) 0;
      c*cos((1.5+1/8)*pi) c*sin((1.5+1/8)*pi) 0;
      c*cos((1.75+1/8)*pi) c*sin((1.75+1/8)*pi) 0; % end of third loop
      ]; 
      
ori.panel0{1}=[1 2 3];
ori.panel0{2}=[1 3 4];
ori.panel0{3}=[1 4 5];
ori.panel0{4}=[1 5 6];
ori.panel0{5}=[1 6 7];
ori.panel0{6}=[1 7 8];
ori.panel0{7}=[1 8 9];
ori.panel0{8}=[1 9 2];

ori.panel0{9}=[2 10 11 3];
ori.panel0{10}=[3 12 13 4];
ori.panel0{11}=[4 14 15 5];
ori.panel0{12}=[5 16 17 6];
ori.panel0{13}=[6 18 19 7];
ori.panel0{14}=[7 20 21 8];
ori.panel0{15}=[8 22 23 9];
ori.panel0{16}=[9 24 25 2];

ori.panel0{17}=[2 26 10];
ori.panel0{18}=[3 11 27];
ori.panel0{19}=[3 27 12];
ori.panel0{20}=[4 13 28];
ori.panel0{21}=[4 28 14];
ori.panel0{22}=[5 15 29];
ori.panel0{23}=[5 29 16];
ori.panel0{24}=[6 17 30];
ori.panel0{25}=[6 30 18];
ori.panel0{26}=[7 19 31];
ori.panel0{27}=[7 31 20];
ori.panel0{28}=[8 21 32];
ori.panel0{29}=[8 32 22];
ori.panel0{30}=[9 23 33];
ori.panel0{31}=[9 33 24];
ori.panel0{32}=[2 25 26];


ori.panel0{33}=[10 34 11];
ori.panel0{34}=[12 35 13];
ori.panel0{35}=[14 36 15];
ori.panel0{36}=[16 37 17];
ori.panel0{37}=[18 38 19];
ori.panel0{38}=[20 39 21];
ori.panel0{39}=[22 40 23];
ori.panel0{40}=[24 41 25];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model
% Define the crease width
ori.creaseWidthMat=zeros(ori.oldCreaseNum,1);

ori.creaseWidthMat(1)=1.5*10^(-3);
ori.creaseWidthMat(2)=1.5*10^(-3);
ori.creaseWidthMat(14)=1.5*10^(-3);
ori.creaseWidthMat(12)=1.5*10^(-3);
ori.creaseWidthMat(10)=1.5*10^(-3);
ori.creaseWidthMat(8)=1.5*10^(-3);
ori.creaseWidthMat(6)=1.5*10^(-3);
ori.creaseWidthMat(4)=1.5*10^(-3);

ori.creaseWidthMat(3)=1.5*10^(-3);
ori.creaseWidthMat(5)=1.5*10^(-3);
ori.creaseWidthMat(7)=1.5*10^(-3);
ori.creaseWidthMat(9)=1.5*10^(-3);
ori.creaseWidthMat(11)=1.5*10^(-3);
ori.creaseWidthMat(13)=1.5*10^(-3);
ori.creaseWidthMat(15)=1.5*10^(-3);
ori.creaseWidthMat(16)=1.5*10^(-3);

ori.creaseWidthMat(41)=1.5*10^(-3);
ori.creaseWidthMat(43)=1.5*10^(-3);
ori.creaseWidthMat(46)=1.5*10^(-3);
ori.creaseWidthMat(49)=1.5*10^(-3);
ori.creaseWidthMat(52)=1.5*10^(-3);
ori.creaseWidthMat(55)=1.5*10^(-3);
ori.creaseWidthMat(58)=1.5*10^(-3);
ori.creaseWidthMat(61)=1.5*10^(-3);

ori.creaseWidthMat(40)=1.5*10^(-3);
ori.creaseWidthMat(17)=1.5*10^(-3);
ori.creaseWidthMat(19)=1.5*10^(-3);
ori.creaseWidthMat(20)=1.5*10^(-3);
ori.creaseWidthMat(22)=1.5*10^(-3);
ori.creaseWidthMat(23)=1.5*10^(-3);
ori.creaseWidthMat(25)=1.5*10^(-3);
ori.creaseWidthMat(26)=1.5*10^(-3);
ori.creaseWidthMat(28)=1.5*10^(-3);
ori.creaseWidthMat(29)=1.5*10^(-3);
ori.creaseWidthMat(31)=1.5*10^(-3);
ori.creaseWidthMat(32)=1.5*10^(-3);
ori.creaseWidthMat(34)=1.5*10^(-3);
ori.creaseWidthMat(35)=1.5*10^(-3);
ori.creaseWidthMat(37)=1.5*10^(-3);
ori.creaseWidthMat(38)=1.5*10^(-3);

ori.creaseWidthMat(18)=1.5*10^(-3);
ori.creaseWidthMat(21)=1.5*10^(-3);
ori.creaseWidthMat(24)=1.5*10^(-3);
ori.creaseWidthMat(27)=1.5*10^(-3);
ori.creaseWidthMat(30)=1.5*10^(-3);
ori.creaseWidthMat(33)=1.5*10^(-3);
ori.creaseWidthMat(36)=1.5*10^(-3);
ori.creaseWidthMat(39)=1.5*10^(-3);

% Compute the meshed geometry
ori.Mesh_CompliantCreaseGeometry()

% Plot the results for inspection
ori.viewAngle1=300;
ori.viewAngle2=30;
ori.displayRange=50*10^(-3); % plotting range
ori.displayRangeRatio=1; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2000*10^6;  
ori.creaseE=20*10^6; 
ori.panelPoisson=0.3; 
ori.creasePoisson=0.3; 
ori.panelThickMat=ones(40,1)*1000*10^(-6);
ori.panelW=1.5*10^-3;

ori.creaseThickMat=zeros(ori.oldCreaseNum,1);
ori.creaseThickMat(1)=300*10^(-6);
ori.creaseThickMat(2)=300*10^(-6);
ori.creaseThickMat(14)=300*10^(-6);
ori.creaseThickMat(12)=300*10^(-6);
ori.creaseThickMat(10)=300*10^(-6);
ori.creaseThickMat(8)=300*10^(-6);
ori.creaseThickMat(6)=300*10^(-6);
ori.creaseThickMat(4)=300*10^(-6);

ori.creaseThickMat(3)=300*10^(-6);
ori.creaseThickMat(5)=300*10^(-6);
ori.creaseThickMat(7)=300*10^(-6);
ori.creaseThickMat(9)=300*10^(-6);
ori.creaseThickMat(11)=300*10^(-6);
ori.creaseThickMat(13)=300*10^(-6);
ori.creaseThickMat(15)=300*10^(-6);
ori.creaseThickMat(16)=300*10^(-6);

ori.creaseThickMat(41)=300*10^(-6);
ori.creaseThickMat(43)=300*10^(-6);
ori.creaseThickMat(46)=300*10^(-6);
ori.creaseThickMat(49)=300*10^(-6);
ori.creaseThickMat(52)=300*10^(-6);
ori.creaseThickMat(55)=300*10^(-6);
ori.creaseThickMat(58)=300*10^(-6);
ori.creaseThickMat(61)=300*10^(-6);

ori.creaseThickMat(40)=300*10^(-6);
ori.creaseThickMat(17)=300*10^(-6);
ori.creaseThickMat(19)=300*10^(-6);
ori.creaseThickMat(20)=300*10^(-6);
ori.creaseThickMat(22)=300*10^(-6);
ori.creaseThickMat(23)=300*10^(-6);
ori.creaseThickMat(25)=300*10^(-6);
ori.creaseThickMat(26)=300*10^(-6);
ori.creaseThickMat(28)=300*10^(-6);
ori.creaseThickMat(29)=300*10^(-6);
ori.creaseThickMat(31)=300*10^(-6);
ori.creaseThickMat(32)=300*10^(-6);
ori.creaseThickMat(34)=300*10^(-6);
ori.creaseThickMat(35)=300*10^(-6);
ori.creaseThickMat(37)=300*10^(-6);
ori.creaseThickMat(38)=300*10^(-6);

ori.creaseThickMat(18)=300*10^(-6);
ori.creaseThickMat(21)=300*10^(-6);
ori.creaseThickMat(24)=300*10^(-6);
ori.creaseThickMat(27)=300*10^(-6);
ori.creaseThickMat(30)=300*10^(-6);
ori.creaseThickMat(33)=300*10^(-6);
ori.creaseThickMat(36)=300*10^(-6);
ori.creaseThickMat(39)=300*10^(-6);

%% Setup the loading controller

selfFold=ControllerSelfFolding;

% Assign zero strain position for creases during self-folding
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);

GlobalRate=0.6;
ratio=-0.4*GlobalRate;
selfFold.targetRotZeroStrain(3)=pi-ratio*pi;
selfFold.targetRotZeroStrain(5)=pi-ratio*pi;
selfFold.targetRotZeroStrain(7)=pi-ratio*pi;
selfFold.targetRotZeroStrain(9)=pi-ratio*pi;
selfFold.targetRotZeroStrain(11)=pi-ratio*pi;
selfFold.targetRotZeroStrain(13)=pi-ratio*pi;
selfFold.targetRotZeroStrain(15)=pi-ratio*pi;
selfFold.targetRotZeroStrain(16)=pi-ratio*pi;

ratio1=-0.4*GlobalRate;
selfFold.targetRotZeroStrain(41)=pi+ratio1*pi;
selfFold.targetRotZeroStrain(43)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(46)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(49)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(52)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(55)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(58)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(61)=pi-ratio1*pi;

selfFold.targetRotZeroStrain(40)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(17)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(19)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(20)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(22)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(23)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(25)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(26)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(28)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(29)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(31)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(32)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(34)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(35)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(37)=pi-ratio1*pi;
selfFold.targetRotZeroStrain(38)=pi-ratio1*pi;

ratio2=-0.4*GlobalRate;
selfFold.targetRotZeroStrain(18)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(21)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(24)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(27)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(30)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(33)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(36)=pi-ratio2*pi;
selfFold.targetRotZeroStrain(39)=pi-ratio2*pi;

selfFold.supp=[1,1,1,1;
               2,1,1,1;
               3,1,1,1;];
  
selfFold.increStep=40;
selfFold.tol=5*10^-5;
selfFold.iterMax=50;

% add a second self fold step
selfFold2=ControllerSelfFolding;
GlobalRate=1;
ratio2=0*GlobalRate;
selfFold2.targetRotZeroStrain=selfFold.targetRotZeroStrain;
selfFold2.targetRotZeroStrain(18)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(21)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(24)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(27)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(30)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(33)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(36)=pi-ratio2*pi;
selfFold2.targetRotZeroStrain(39)=pi-ratio2*pi;

selfFold2.supp=[1,1,1,1;
               2,1,1,1;
               3,1,1,1;];
  
selfFold2.increStep=40;
selfFold2.tol=5*10^-5;
selfFold2.iterMax=50;

% Here shows a way to stack different classes in the solver to achieve
% multi-step loading of the origami structure. 
ori.loadingController{1}={"SelfFold",selfFold};
ori.loadingController{2}={"SelfFold",selfFold2};


%% Solving the model
ori.Solver_Solve();

