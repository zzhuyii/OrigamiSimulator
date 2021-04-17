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

% thickness of gold and su-8
tg=0.2*10^-6;
ts=0.9*10^-6;
% density of SU-8
rhoSU8=1200;
% thickness of panel
tpanel=21*10^-6;
% length of panel
panelL=450*10^(-6);
% width of creases
W=200*10^(-6);
% undercut of XeF2 etch
undercut=60*10^-6;
% residual fold during fabrication
residualFold=-15;
% power input (mW)
qload=10;


ori.node0=[0 0 0;
      panelL 0 0;
      2*panelL+W 0 0;
      3*panelL+1.5*W 0 0;
      0 panelL 0;
      panelL panelL 0;
      2*panelL+W panelL 0;
      3*panelL+1.5*W panelL 0;
      0*panelL,-0.6*panelL,-undercut;
      0*panelL,1.6*panelL,-undercut;
      4*panelL,-0.6*panelL,-undercut;
      4*panelL,1.6*panelL,-undercut;];
  
ori.panel0{1}=[1 2 6 5];
ori.panel0{2}=[2 3 7 6];
ori.panel0{3}=[3 4 8 7];
ori.panel0{4}=[9 11 12 10];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;
ori.creaseWidthVec(6)=W;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=2*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;



%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[500*10^(-6);tpanel;tpanel;500*10^(-6)]; 
ori.panelW=W;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=1000;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(tg+ts);
ori.creaseThickVec(6)=(tg+ts);


%% setup panel contact information

ori.contactOpen=1;
ori.ke=0.0001;
ori.d0edge=20*(10^(-6));
ori.d0center=20*(10^(-6));


%% Assign Thermal Properties

ori.panelThermalConductVec = [1.3;0.3;0.3;1.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=600*10^(-6); 



%% Setup the loading controller

% applying the residual stress
selfFold=ControllerSelfFolding;

% Assign zero strain position for creases during self-folding
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);

ratio=residualFold/180;
selfFold.targetRotZeroStrain(3)=pi+ratio*pi;
selfFold.targetRotZeroStrain(6)=pi+ratio*pi;

selfFold.supp=[1,0,0,1;
      2,0,0,1;
      3,0,1,1;
      4,1,1,1;
      13,1,1,1;
      14,1,1,1;
      15,1,1,1;
      16,1,1,1;];
selfFold.increStep=50;
selfFold.tol=4*10^-6; 
selfFold.iterMax=25;
selfFold.videoOpen=0;


% applying the gravity loading
nr=ControllerNRLoading;

nr.increStep=80;
nr.tol=2*10^-6;
nr.iterMax=50;

nr.supp=[1,0,0,1;
      2,0,0,1;
      3,0,1,1;
      4,1,1,1;
      13,1,1,1;
      14,1,1,1;
      15,1,1,1;
      16,1,1,1;];

loadMag=0.9*(panelL*panelL*tpanel*rhoSU8*9.8)/nr.increStep/6;
nr.load=[5,0,0,-loadMag;
      6,0,0,-loadMag;
      7,0,0,-loadMag;
      8,0,0,-loadMag;
      24,0,0,-2*loadMag;
      9,0,0,-loadMag;
      10,0,0,-loadMag;
      11,0,0,-loadMag;
      12,0,0,-loadMag;
      25,0,0,-2*loadMag;];
nr.videoOpen=0;

ori.loadingController{1}={"SelfFold",selfFold};
ori.loadingController{2}={"NR",nr};
ori.Solver_Solve();

ori.continuingLoading=1;
gravityShape=ori.newNode+ori.currentU;

% applying the thermal loading
thermal=ControllerThermalLoading;

thermal.thermalStep=250;
thermal.tol=1*10^-6; 

thermal.supp=[1,0,0,1;
      2,0,0,1;
      3,0,1,1;
      4,1,1,1;
      13,1,1,1;
      14,1,1,1;
      15,1,1,1;
      16,1,1,1;];

thermal.thermalBoundaryPanelVec=[4];
thermal.roomTempNode=[1;4];

thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
thermal.deltaAlpha(3)=(52-14)*10^(-6); 
thermal.deltaAlpha(6)=(52-14)*10^(-6); 

thermal.Emat1=0.7*79*10^9; 
thermal.Emat2=2*10^9;
thermal.tmat1=tg;
thermal.tmat2=ts;
thermal.videoOpen=0; % close the animation

% the target loading of crease heating
thermal.targetCreaseHeating=[3,qload/1000;6,qload/1000];

ori.loadingController={};
ori.loadingController{1}={"ThermalLoading",thermal};

ori.Solver_Solve();


%% Solving the folding

fold1His=zeros(thermal.thermalStep,1);
fold2His=zeros(thermal.thermalStep,1);

for i=1:thermal.thermalStep
    
    node1=squeeze(gravityShape(1,:))'+squeeze(thermal.Uhis(i,1,:));
    node2=squeeze(gravityShape(2,:))'+squeeze(thermal.Uhis(i,2,:));
    node3=squeeze(gravityShape(5,:))'+squeeze(thermal.Uhis(i,5,:));
    node4=squeeze(gravityShape(6,:))'+squeeze(thermal.Uhis(i,6,:));
    vec1=node2-node1;
    vec2=node4-node3;
    rotation=dot(vec1,vec2)/norm(vec1)/norm(vec2);
    rotation=sign(node4(3)-node1(3))*acos(rotation);
    fold1His(i)=rotation;
    
    node1=squeeze(gravityShape(5,:))'+squeeze(thermal.Uhis(i,5,:));
    node2=squeeze(gravityShape(6,:))'+squeeze(thermal.Uhis(i,6,:));
    node3=squeeze(gravityShape(9,:))'+squeeze(thermal.Uhis(i,9,:));
    node4=squeeze(gravityShape(10,:))'+squeeze(thermal.Uhis(i,10,:));
    vec1=node2-node1;
    vec2=node4-node3;
    rotation=dot(vec1,vec2)/norm(vec1)/norm(vec2);
    fold2His(i)=sign(node4(3)-(node2(3)-(node1(3)-node2(3))))*acos(rotation);    
    
end


for i=2:thermal.thermalStep

    fold1His(i)=fold1His(i)-fold1His(1);
    fold2His(i)=fold2His(i)-fold2His(1);
    
end
fold1His(1)=0;
fold2His(1)=0;