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
a=15*10^-3;
b=3*10^-3;
l=30*10^-3;

bw=1.5*10^-3;
bm=1.5*10^-3;

yshift=0*a;
syshift=0*a;

minorshift=0.2*l;
wshift=0.5*l;
oshift=2*l;
mshift=6.2*l;
pshift=6.7*l;
sshift=10.5*l;



% thickness of two layer
t1=0.02*10^-3;
t2=0.02*10^-3;

% thickness of panel
tpanel=1*10^-3;

% width of crease
W=2*10^-3;

% power input (W)
qload=0.1;
qloadw=0.2;
qloado=0.10;
qloadm=0.2;
qloadp=0.10;


ori.node0=[minorshift -a 0;
           minorshift 0 0;
           minorshift a 0;
           minorshift+l-b -a 0;
           minorshift+l+b 0 0;
           minorshift+l-b a 0;
           minorshift+2*l+b -a 0;
           minorshift+2*l-b 0 0;
           minorshift+2*l+b a 0;
           minorshift+3*l+b -a 0;
           minorshift+3*l-b 0 0;
           minorshift+3*l+b a 0;
           minorshift+4*l-b -a 0;
           minorshift+4*l+b 0 0;
           minorshift+4*l-b a 0;
           minorshift+5*l -a 0;
           minorshift+5*l 0 0;
           minorshift+5*l a 0;
           
           wshift yshift+-a 0;
           wshift yshift+0 0;
           wshift yshift+a 0;
           wshift+2*l+bw yshift+-a 0;
           wshift+2*l-bw yshift+0 0;
           wshift+2*l+bw yshift+a 0;
           wshift+3*l+bw yshift+-a 0;
           wshift+3*l-bw yshift+0 0;
           wshift+3*l+bw yshift+a 0;
           wshift+4*l+bw yshift+-a 0;
           wshift+4*l-bw yshift+0 0;
           wshift+4*l+bw yshift+a 0;
           wshift+6*l yshift+-a 0;
           wshift+6*l yshift+0 0;
           wshift+6*l yshift+a 0;
           
           
           oshift -a 0;
           oshift 0 0;
           oshift a 0;
           oshift+0.5*l+b -a 0;
           oshift+0.5*l-b 0 0;
           oshift+0.5*l+b a 0;
           oshift+2.5*l-b -a 0;
           oshift+2.5*l+b 0 0;
           oshift+2.5*l-b a 0;
           oshift+3.5*l+b -a 0;
           oshift+3.5*l-b 0 0;
           oshift+3.5*l+b a 0;
           oshift+5.5*l-b -a 0;
           oshift+5.5*l+b 0 0;
           oshift+5.5*l-b a 0;
           oshift+6*l -a 0;
           oshift+6*l 0 0;
           oshift+6*l a 0;
                      
           mshift yshift+-a 0;
           mshift yshift+0 0;
           mshift yshift+a 0;
           mshift+2*l+bm yshift+-a 0;
           mshift+2*l-bm yshift+0 0;
           mshift+2*l+bm yshift+a 0;
           mshift+3*l+bm yshift+-a 0;
           mshift+3*l-bm yshift+0 0;
           mshift+3*l+bm yshift+a 0;
           mshift+4*l+bm yshift+-a 0;
           mshift+4*l-bm yshift+0 0;
           mshift+4*l+bm yshift+a 0;
           mshift+6*l yshift+-a 0;
           mshift+6*l yshift+0 0;
           mshift+6*l yshift+a 0;
           
           pshift -a 2*l;
           pshift 0 2*l;
           pshift a 2*l;
           pshift+2*l+b -a 2*l;
           pshift+2*l-b 0 2*l;
           pshift+2*l+b a 2*l;
           pshift+3*l-b -a 2*l;
           pshift+3*l+b 0 2*l;
           pshift+3*l-b a 2*l;
           pshift+4*l+b -a 2*l;
           pshift+4*l-b 0 2*l;
           pshift+4*l+b a 2*l;
           pshift+4.7*l -a 2*l;
           pshift+4.7*l 0 2*l;
           pshift+4.7*l a 2*l;
           
           sshift syshift+-a 0;
           sshift syshift+0 0;
           sshift syshift+a 0;
           sshift+l-b syshift+-a 0;
           sshift+l+b syshift+0 0;
           sshift+l-b syshift+a 0;
           sshift+2*l+b syshift+-a 0;
           sshift+2*l-b syshift+0 0;
           sshift+2*l+b syshift+a 0;
           sshift+3*l+b syshift+-a 0;
           sshift+3*l-b syshift+0 0;
           sshift+3*l+b syshift+a 0;
           sshift+4*l-b syshift+-a 0;
           sshift+4*l+b syshift+0 0;
           sshift+4*l-b syshift+a 0;
           sshift+5*l syshift+-a 0;
           sshift+5*l syshift+0 0;
           sshift+5*l syshift+a 0;
           ];
       
ori.node0(:,1)=ori.node0(:,1)-3*l;
  
ori.panel0{1}=[1 4 5 2];
ori.panel0{2}=[2 5 6 3];
ori.panel0{3}=[4 7 8 5];
ori.panel0{4}=[5 8 9 6];
ori.panel0{5}=[7 10 11 8];
ori.panel0{6}=[8 11 12 9];
ori.panel0{7}=[10 13 14 11];
ori.panel0{8}=[11 14 15 12];
ori.panel0{9}=[13 16 17 14];
ori.panel0{10}=[14 17 18 15];

ori.panel0{11}=18+[1 4 5 2];
ori.panel0{12}=18+[2 5 6 3];
ori.panel0{13}=18+[4 7 8 5];
ori.panel0{14}=18+[5 8 9 6];
ori.panel0{15}=18+[7 10 11 8];
ori.panel0{16}=18+[8 11 12 9];
ori.panel0{17}=18+[10 13 14 11];
ori.panel0{18}=18+[11 14 15 12];

ori.panel0{19}=33+[1 4 5 2];
ori.panel0{20}=33+[2 5 6 3];
ori.panel0{21}=33+[4 7 8 5];
ori.panel0{22}=33+[5 8 9 6];
ori.panel0{23}=33+[7 10 11 8];
ori.panel0{24}=33+[8 11 12 9];
ori.panel0{25}=33+[10 13 14 11];
ori.panel0{26}=33+[11 14 15 12];
ori.panel0{27}=33+[13 16 17 14];
ori.panel0{28}=33+[14 17 18 15];

ori.panel0{29}=51+[1 4 5 2];
ori.panel0{30}=51+[2 5 6 3];
ori.panel0{31}=51+[4 7 8 5];
ori.panel0{32}=51+[5 8 9 6];
ori.panel0{33}=51+[7 10 11 8];
ori.panel0{34}=51+[8 11 12 9];
ori.panel0{35}=51+[10 13 14 11];
ori.panel0{36}=51+[11 14 15 12];

ori.panel0{37}=66+[1 4 5 2];
ori.panel0{38}=66+[2 5 6 3];
ori.panel0{39}=66+[4 7 8 5];
ori.panel0{40}=66+[5 8 9 6];
ori.panel0{41}=66+[7 10 11 8];
ori.panel0{42}=66+[8 11 12 9];
ori.panel0{43}=66+[10 13 14 11];
ori.panel0{44}=66+[11 14 15 12];

ori.panel0{37}=66+[1 4 5 2];
ori.panel0{38}=66+[2 5 6 3];
ori.panel0{39}=66+[4 7 8 5];
ori.panel0{40}=66+[5 8 9 6];
ori.panel0{41}=66+[7 10 11 8];
ori.panel0{42}=66+[8 11 12 9];
ori.panel0{43}=66+[10 13 14 11];
ori.panel0{44}=66+[11 14 15 12];

ori.panel0{45}=81+[1 4 5 2];
ori.panel0{46}=81+[2 5 6 3];
ori.panel0{47}=81+[4 7 8 5];
ori.panel0{48}=81+[5 8 9 6];
ori.panel0{49}=81+[7 10 11 8];
ori.panel0{50}=81+[8 11 12 9];
ori.panel0{51}=81+[10 13 14 11];
ori.panel0{52}=81+[11 14 15 12];
ori.panel0{53}=81+[13 16 17 14];
ori.panel0{54}=81+[14 17 18 15];


numPanel=54;

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();

% Plot the results for inspection
ori.viewAngle1=-7;
ori.viewAngle2=4;
ori.displayRange=400*10^(-3); % plotting range
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
ori.creaseWidthVec(19)=W;
ori.creaseWidthVec(21)=W;
ori.creaseWidthVec(25)=W;

ori.creaseWidthVec(31)=W;
ori.creaseWidthVec(30)=W;
ori.creaseWidthVec(33)=W;
ori.creaseWidthVec(37)=W;
ori.creaseWidthVec(36)=W;
ori.creaseWidthVec(38)=W;
ori.creaseWidthVec(42)=W;
ori.creaseWidthVec(41)=W;
ori.creaseWidthVec(43)=W;
ori.creaseWidthVec(47)=W;

ori.creaseWidthVec(49+3)=W;
ori.creaseWidthVec(49+4)=W;
ori.creaseWidthVec(49+6)=W;
ori.creaseWidthVec(49+10)=W;
ori.creaseWidthVec(49+9)=W;
ori.creaseWidthVec(49+11)=W;
ori.creaseWidthVec(49+15)=W;
ori.creaseWidthVec(49+14)=W;
ori.creaseWidthVec(49+16)=W;
ori.creaseWidthVec(49+20)=W;
ori.creaseWidthVec(49+19)=W;
ori.creaseWidthVec(49+21)=W;
ori.creaseWidthVec(49+25)=W;

ori.creaseWidthVec(76+3)=W;
ori.creaseWidthVec(76+4)=W;
ori.creaseWidthVec(76+6)=W;
ori.creaseWidthVec(76+10)=W;
ori.creaseWidthVec(76+9)=W;
ori.creaseWidthVec(76+11)=W;
ori.creaseWidthVec(76+15)=W;
ori.creaseWidthVec(76+14)=W;
ori.creaseWidthVec(76+16)=W;
ori.creaseWidthVec(76+20)=W;

ori.creaseWidthVec(98+3)=W;
ori.creaseWidthVec(98+4)=W;
ori.creaseWidthVec(98+6)=W;
ori.creaseWidthVec(98+10)=W;
ori.creaseWidthVec(98+9)=W;
ori.creaseWidthVec(98+11)=W;
ori.creaseWidthVec(98+15)=W;
ori.creaseWidthVec(98+14)=W;
ori.creaseWidthVec(98+16)=W;
ori.creaseWidthVec(98+20)=W;

ori.creaseWidthVec(120+3)=W;
ori.creaseWidthVec(120+4)=W;
ori.creaseWidthVec(120+6)=W;
ori.creaseWidthVec(120+10)=W;
ori.creaseWidthVec(120+9)=W;
ori.creaseWidthVec(120+11)=W;
ori.creaseWidthVec(120+15)=W;
ori.creaseWidthVec(120+14)=W;
ori.creaseWidthVec(120+16)=W;
ori.creaseWidthVec(120+20)=W;
ori.creaseWidthVec(120+19)=W;
ori.creaseWidthVec(120+21)=W;
ori.creaseWidthVec(120+25)=W;

% Compute the meshed geometry
ori.Mesh_Mesh()




ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

ori.panelThickVec=tpanel*ones(numPanel,1); 
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
ori.creaseThickVec(19)=t1+t2;
ori.creaseThickVec(21)=t1+t2;
ori.creaseThickVec(25)=t1+t2;

ori.creaseThickVec(31)=t1+t2;
ori.creaseThickVec(30)=t1+t2;
ori.creaseThickVec(33)=t1+t2;
ori.creaseThickVec(37)=t1+t2;
ori.creaseThickVec(38)=t1+t2;
ori.creaseThickVec(36)=t1+t2;
ori.creaseThickVec(42)=t1+t2;
ori.creaseThickVec(41)=t1+t2;
ori.creaseThickVec(43)=t1+t2;
ori.creaseThickVec(47)=t1+t2;

ori.creaseThickVec(49+3)=t1+t2;
ori.creaseThickVec(49+4)=t1+t2;
ori.creaseThickVec(49+6)=t1+t2;
ori.creaseThickVec(49+10)=t1+t2;
ori.creaseThickVec(49+9)=t1+t2;
ori.creaseThickVec(49+11)=t1+t2;
ori.creaseThickVec(49+15)=t1+t2;
ori.creaseThickVec(49+14)=t1+t2;
ori.creaseThickVec(49+16)=t1+t2;
ori.creaseThickVec(49+20)=t1+t2;
ori.creaseThickVec(49+19)=t1+t2;
ori.creaseThickVec(49+21)=t1+t2;
ori.creaseThickVec(49+25)=t1+t2;

ori.creaseThickVec(76+3)=t1+t2;
ori.creaseThickVec(76+4)=t1+t2;
ori.creaseThickVec(76+6)=t1+t2;
ori.creaseThickVec(76+10)=t1+t2;
ori.creaseThickVec(76+9)=t1+t2;
ori.creaseThickVec(76+11)=t1+t2;
ori.creaseThickVec(76+15)=t1+t2;
ori.creaseThickVec(76+14)=t1+t2;
ori.creaseThickVec(76+16)=t1+t2;
ori.creaseThickVec(76+20)=t1+t2;

ori.creaseThickVec(98+3)=t1+t2;
ori.creaseThickVec(98+4)=t1+t2;
ori.creaseThickVec(98+6)=t1+t2;
ori.creaseThickVec(98+10)=t1+t2;
ori.creaseThickVec(98+9)=t1+t2;
ori.creaseThickVec(98+11)=t1+t2;
ori.creaseThickVec(98+15)=t1+t2;
ori.creaseThickVec(98+14)=t1+t2;
ori.creaseThickVec(98+16)=t1+t2;
ori.creaseThickVec(98+20)=t1+t2;

ori.creaseThickVec(120+3)=t1+t2;
ori.creaseThickVec(120+4)=t1+t2;
ori.creaseThickVec(120+6)=t1+t2;
ori.creaseThickVec(120+10)=t1+t2;
ori.creaseThickVec(120+9)=t1+t2;
ori.creaseThickVec(120+11)=t1+t2;
ori.creaseThickVec(120+15)=t1+t2;
ori.creaseThickVec(120+14)=t1+t2;
ori.creaseThickVec(120+16)=t1+t2;
ori.creaseThickVec(120+20)=t1+t2;
ori.creaseThickVec(120+19)=t1+t2;
ori.creaseThickVec(120+21)=t1+t2;
ori.creaseThickVec(120+25)=t1+t2;


%% Assign Thermal Properties

ori.panelThermalConductVec = 0.3*ones(numPanel,1); 

ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=0.02; 



%% Setup the loading controller

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=50;
% I am using 250 when comparing the efficiency

thermal.tol=5*10^-5; 

thermal.supp=[1,1,1,1;
      2,1,1,1;
      7,1,0,1;
      8,1,0,1;
      
      49,1,1,1;
      58,0,1,1;
      56,1,0,1;
      63,0,0,1;
      
      89,1,1,1;
      90,1,1,1;
      96,1,0,1;
      95,1,0,1;
      
      112+1,1,1,1;
      112+26,0,1,1;
      112+8,1,0,1;
      112+31,0,0,1;
      
      153,1,1,1;
      154,0,1,1;
      160,1,0,1;
      159,0,0,1;
      
      177,1,1,1;
      178,1,1,1;
      183,1,0,1;
      184,1,0,1;
      ];

thermal.thermalBoundaryPanelVec=[];
thermal.roomTempNode=[];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
deltaAlpha=1000*10^(-6);
thermal.deltaAlpha(3)= deltaAlpha;
thermal.deltaAlpha(4)= 0;
thermal.deltaAlpha(6)= deltaAlpha;
thermal.deltaAlpha(10)= 0;
thermal.deltaAlpha(9)= deltaAlpha;
thermal.deltaAlpha(11)= deltaAlpha;
thermal.deltaAlpha(15)= 0;
thermal.deltaAlpha(14)= -deltaAlpha;
thermal.deltaAlpha(16)= -deltaAlpha;
thermal.deltaAlpha(20)= 0;
thermal.deltaAlpha(19)= -deltaAlpha;
thermal.deltaAlpha(21)= -deltaAlpha;
thermal.deltaAlpha(25)= 0;

thermal.deltaAlpha(31)= deltaAlpha;
thermal.deltaAlpha(30)= 0;
thermal.deltaAlpha(33)= deltaAlpha;
thermal.deltaAlpha(37)= 0;
thermal.deltaAlpha(36)= -deltaAlpha;
thermal.deltaAlpha(38)= -deltaAlpha;
thermal.deltaAlpha(42)= 0;
thermal.deltaAlpha(41)= deltaAlpha;
thermal.deltaAlpha(43)= deltaAlpha;
thermal.deltaAlpha(43)= 0;

thermal.deltaAlpha(49+3)= deltaAlpha;
thermal.deltaAlpha(49+4)= 0;
thermal.deltaAlpha(49+6)= deltaAlpha;
thermal.deltaAlpha(49+10)= 0;
thermal.deltaAlpha(49+9)= deltaAlpha;
thermal.deltaAlpha(49+11)= deltaAlpha;
thermal.deltaAlpha(49+15)= 0;
thermal.deltaAlpha(49+14)= deltaAlpha;
thermal.deltaAlpha(49+16)= deltaAlpha;
thermal.deltaAlpha(49+20)= 0;
thermal.deltaAlpha(49+19)= deltaAlpha;
thermal.deltaAlpha(49+21)= deltaAlpha;
thermal.deltaAlpha(49+25)= 0;

thermal.deltaAlpha(76+3)= -deltaAlpha;
thermal.deltaAlpha(76+4)= 0;
thermal.deltaAlpha(76+6)= -deltaAlpha;
thermal.deltaAlpha(76+10)= 0;
thermal.deltaAlpha(76+9)= deltaAlpha;
thermal.deltaAlpha(76+11)= deltaAlpha;
thermal.deltaAlpha(76+15)= 0;
thermal.deltaAlpha(76+14)= -deltaAlpha;
thermal.deltaAlpha(76+16)= -deltaAlpha;
thermal.deltaAlpha(76+20)= 0;

thermal.deltaAlpha(98+3)= -deltaAlpha;
thermal.deltaAlpha(98+4)= 0;
thermal.deltaAlpha(98+6)= -deltaAlpha;
thermal.deltaAlpha(98+10)= 0;
thermal.deltaAlpha(98+9)= -deltaAlpha;
thermal.deltaAlpha(98+11)= -deltaAlpha;
thermal.deltaAlpha(98+15)= 0;
thermal.deltaAlpha(98+14)= -deltaAlpha;
thermal.deltaAlpha(98+16)= -deltaAlpha;
thermal.deltaAlpha(98+20)= 0;

thermal.deltaAlpha(120+3)= deltaAlpha;
thermal.deltaAlpha(120+4)= 0;
thermal.deltaAlpha(120+6)= deltaAlpha;
thermal.deltaAlpha(120+10)= 0;
thermal.deltaAlpha(120+9)= deltaAlpha;
thermal.deltaAlpha(120+11)= deltaAlpha;
thermal.deltaAlpha(120+15)= 0;
thermal.deltaAlpha(120+14)= -deltaAlpha;
thermal.deltaAlpha(120+16)= -deltaAlpha;
thermal.deltaAlpha(120+20)= 0;
thermal.deltaAlpha(120+19)= -deltaAlpha;
thermal.deltaAlpha(120+21)= -deltaAlpha;
thermal.deltaAlpha(120+25)= 0;


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
    20,qload;
    19,qload;
    21,qload;
    25,qload;
    
    31,qloadw;
    30,qloadw;
    33,qloadw;
    37,qloadw;
    36,qloadw;    
    38,qloadw;
    42,qloadw;
    41,qloadw;
    43,qloadw;
    47,qloadw;
    
    49+3,qloado;
    49+4,qloado;
    49+6,qloado;
    49+10,qloado;
    49+9,qloado;
    49+11,qloado;
    49+15,qloado;
    49+14,qloado;    
    49+16,qloado;
    49+20,qloado;
    49+19,qloado;
    49+21,qloado;
    49+25,qloado;
    
    76+3,qloadm;
    76+4,qloadm;
    76+6,qloadm;
    76+10,qloadm;
    76+9,qloadm;
    76+11,qloadm;
    76+15,qloadm;
    76+14,qloadm;    
    76+16,qloadm;
    76+20,qloadm;
    
    98+3,qloadp;
    98+4,qloadp;
    98+6,qloadp;
    98+10,qloadp;
    98+9,qloadp;
    98+11,qloadp;
    98+15,qloadp;
    98+14,qloadp;    
    98+16,qloadp;
    98+20,qloadp;
    
    
    120+3,qload;
    120+4,qload;
    120+6,qload;
    120+10,qload;
    120+9,qload;
    120+11,qload;
    120+15,qload;
    120+14,qload;    
    120+16,qload;
    120+20,qload;
    120+19,qload;
    120+21,qload;
    120+25,qload;
    ];

ori.loadingController{1}={"ElectroThermal",thermal};
ori.Solver_Solve();

toc
