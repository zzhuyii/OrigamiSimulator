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
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);

ori.creaseWidthVec(1)=1.5*10^(-3);
ori.creaseWidthVec(2)=1.5*10^(-3);
ori.creaseWidthVec(14)=1.5*10^(-3);
ori.creaseWidthVec(12)=1.5*10^(-3);
ori.creaseWidthVec(10)=1.5*10^(-3);
ori.creaseWidthVec(8)=1.5*10^(-3);
ori.creaseWidthVec(6)=1.5*10^(-3);
ori.creaseWidthVec(4)=1.5*10^(-3);

ori.creaseWidthVec(3)=1.5*10^(-3);
ori.creaseWidthVec(5)=1.5*10^(-3);
ori.creaseWidthVec(7)=1.5*10^(-3);
ori.creaseWidthVec(9)=1.5*10^(-3);
ori.creaseWidthVec(11)=1.5*10^(-3);
ori.creaseWidthVec(13)=1.5*10^(-3);
ori.creaseWidthVec(15)=1.5*10^(-3);
ori.creaseWidthVec(16)=1.5*10^(-3);

ori.creaseWidthVec(41)=1.5*10^(-3);
ori.creaseWidthVec(43)=1.5*10^(-3);
ori.creaseWidthVec(46)=1.5*10^(-3);
ori.creaseWidthVec(49)=1.5*10^(-3);
ori.creaseWidthVec(52)=1.5*10^(-3);
ori.creaseWidthVec(55)=1.5*10^(-3);
ori.creaseWidthVec(58)=1.5*10^(-3);
ori.creaseWidthVec(61)=1.5*10^(-3);

ori.creaseWidthVec(40)=1.5*10^(-3);
ori.creaseWidthVec(17)=1.5*10^(-3);
ori.creaseWidthVec(19)=1.5*10^(-3);
ori.creaseWidthVec(20)=1.5*10^(-3);
ori.creaseWidthVec(22)=1.5*10^(-3);
ori.creaseWidthVec(23)=1.5*10^(-3);
ori.creaseWidthVec(25)=1.5*10^(-3);
ori.creaseWidthVec(26)=1.5*10^(-3);
ori.creaseWidthVec(28)=1.5*10^(-3);
ori.creaseWidthVec(29)=1.5*10^(-3);
ori.creaseWidthVec(31)=1.5*10^(-3);
ori.creaseWidthVec(32)=1.5*10^(-3);
ori.creaseWidthVec(34)=1.5*10^(-3);
ori.creaseWidthVec(35)=1.5*10^(-3);
ori.creaseWidthVec(37)=1.5*10^(-3);
ori.creaseWidthVec(38)=1.5*10^(-3);

ori.creaseWidthVec(18)=1.5*10^(-3);
ori.creaseWidthVec(21)=1.5*10^(-3);
ori.creaseWidthVec(24)=1.5*10^(-3);
ori.creaseWidthVec(27)=1.5*10^(-3);
ori.creaseWidthVec(30)=1.5*10^(-3);
ori.creaseWidthVec(33)=1.5*10^(-3);
ori.creaseWidthVec(36)=1.5*10^(-3);
ori.creaseWidthVec(39)=1.5*10^(-3);

% Compute the meshed geometry
ori.Mesh_Mesh()

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
ori.panelThickVec=ones(40,1)*3*10^(-3);
ori.panelW=1.5*10^-3;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(1)=0.5*10^(-3);
ori.creaseThickVec(2)=0.5*10^(-3);
ori.creaseThickVec(14)=0.5*10^(-3);
ori.creaseThickVec(12)=0.5*10^(-3);
ori.creaseThickVec(10)=0.5*10^(-3);
ori.creaseThickVec(8)=0.5*10^(-3);
ori.creaseThickVec(6)=0.5*10^(-3);
ori.creaseThickVec(4)=0.5*10^(-3);

ori.creaseThickVec(3)=0.5*10^(-3);
ori.creaseThickVec(5)=0.5*10^(-3);
ori.creaseThickVec(7)=0.5*10^(-3);
ori.creaseThickVec(9)=0.5*10^(-3);
ori.creaseThickVec(11)=0.5*10^(-3);
ori.creaseThickVec(13)=0.5*10^(-3);
ori.creaseThickVec(15)=0.5*10^(-3);
ori.creaseThickVec(16)=0.5*10^(-3);

ori.creaseThickVec(41)=0.5*10^(-3);
ori.creaseThickVec(43)=0.5*10^(-3);
ori.creaseThickVec(46)=0.5*10^(-3);
ori.creaseThickVec(49)=0.5*10^(-3);
ori.creaseThickVec(52)=0.5*10^(-3);
ori.creaseThickVec(55)=0.5*10^(-3);
ori.creaseThickVec(58)=0.5*10^(-3);
ori.creaseThickVec(61)=0.5*10^(-3);

ori.creaseThickVec(40)=0.5*10^(-3);
ori.creaseThickVec(17)=0.5*10^(-3);
ori.creaseThickVec(19)=0.5*10^(-3);
ori.creaseThickVec(20)=0.5*10^(-3);
ori.creaseThickVec(22)=0.5*10^(-3);
ori.creaseThickVec(23)=0.5*10^(-3);
ori.creaseThickVec(25)=0.5*10^(-3);
ori.creaseThickVec(26)=0.5*10^(-3);
ori.creaseThickVec(28)=0.5*10^(-3);
ori.creaseThickVec(29)=0.5*10^(-3);
ori.creaseThickVec(31)=0.5*10^(-3);
ori.creaseThickVec(32)=0.5*10^(-3);
ori.creaseThickVec(34)=0.5*10^(-3);
ori.creaseThickVec(35)=0.5*10^(-3);
ori.creaseThickVec(37)=0.5*10^(-3);
ori.creaseThickVec(38)=0.5*10^(-3);

ori.creaseThickVec(18)=0.5*10^(-3);
ori.creaseThickVec(21)=0.5*10^(-3);
ori.creaseThickVec(24)=0.5*10^(-3);
ori.creaseThickVec(27)=0.5*10^(-3);
ori.creaseThickVec(30)=0.5*10^(-3);
ori.creaseThickVec(33)=0.5*10^(-3);
ori.creaseThickVec(36)=0.5*10^(-3);
ori.creaseThickVec(39)=0.5*10^(-3);


%% Assign Thermal Properties

ori.panelThermalConductVec = ones(40,1)*1.3; 
ori.creaseThermalConduct=1.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=15*10^(-3); 


%% Setup the loading controller

% applying the thermal loading
thermal=ControllerElectroThermalFolding;

thermal.thermalStep=40;
thermal.tol=8*10^-4; 

thermal.supp=[1,1,1,1;
               2,1,1,1;
               3,1,1,1;];

thermal.thermalBoundaryPanelVec=[];
thermal.roomTempNode=[];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);

deltaAlpha=10^-3;

thermal.deltaAlpha(3)=deltaAlpha;
thermal.deltaAlpha(5)=deltaAlpha;
thermal.deltaAlpha(7)=deltaAlpha;
thermal.deltaAlpha(9)=deltaAlpha;
thermal.deltaAlpha(11)=deltaAlpha; 
thermal.deltaAlpha(13)=deltaAlpha; 
thermal.deltaAlpha(15)=deltaAlpha; 
thermal.deltaAlpha(16)=deltaAlpha;

thermal.deltaAlpha(41)=-deltaAlpha;
thermal.deltaAlpha(43)=-deltaAlpha;
thermal.deltaAlpha(46)=-deltaAlpha;
thermal.deltaAlpha(49)=-deltaAlpha;
thermal.deltaAlpha(52)=-deltaAlpha;
thermal.deltaAlpha(55)=-deltaAlpha;
thermal.deltaAlpha(58)=-deltaAlpha;
thermal.deltaAlpha(61)=-deltaAlpha;

thermal.deltaAlpha(40)=deltaAlpha;
thermal.deltaAlpha(17)=deltaAlpha;
thermal.deltaAlpha(19)=deltaAlpha;
thermal.deltaAlpha(20)=deltaAlpha;
thermal.deltaAlpha(22)=deltaAlpha;
thermal.deltaAlpha(23)=deltaAlpha;
thermal.deltaAlpha(25)=deltaAlpha;
thermal.deltaAlpha(26)=deltaAlpha;
thermal.deltaAlpha(28)=deltaAlpha;
thermal.deltaAlpha(29)=deltaAlpha;
thermal.deltaAlpha(31)=deltaAlpha;
thermal.deltaAlpha(32)=deltaAlpha;
thermal.deltaAlpha(34)=deltaAlpha;
thermal.deltaAlpha(35)=deltaAlpha;
thermal.deltaAlpha(37)=deltaAlpha;
thermal.deltaAlpha(38)=deltaAlpha;

thermal.deltaAlpha(18)=deltaAlpha;
thermal.deltaAlpha(21)=deltaAlpha;
thermal.deltaAlpha(24)=deltaAlpha;
thermal.deltaAlpha(27)=deltaAlpha;
thermal.deltaAlpha(30)=deltaAlpha;
thermal.deltaAlpha(33)=deltaAlpha;
thermal.deltaAlpha(36)=deltaAlpha;
thermal.deltaAlpha(39)=deltaAlpha;


thermal.Emat1=2*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=0.25*10^-3;
thermal.tmat2=0.25*10^-3;
thermal.videoOpen=0;

% the target loading of crease heating
qload=0.4;
thermal.targetCreaseHeating=[
    3,qload;
    5,qload;
    7,qload;
    11,qload;
    13,qload;
    15,qload;
    16,qload;
    41,qload;
    43,qload;
    46,qload;
    49,qload;
    52,qload;
    55,qload;
    58,qload;
    61,qload;
    40,qload;
    17,qload;
    19,qload;
    20,qload;
    22,qload;
    23,qload;
    25,qload;
    26,qload;
    28,qload;
    29,qload;
    31,qload;
    32,qload;
    34,qload;
    35,qload;
    37,qload;
    38,qload;
    18,qload;
    21,qload;
    24,qload;
    27,qload;
    30,qload;
    33,qload;
    36,qload;
    39,qload;];

ori.loadingController{1}={"ElectroThermal",thermal};


%% Solving the model
ori.Solver_Solve();
toc
