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

% Geometry control
barArea=0.000001;
hingeStiff=100;
centerSeparation=0.4;

% Set up Solver
ori=OrigamiSolver;
ori.compliantCreaseOpen=0;
ori.continuingLoading=1;
ori.panelInnerBarStart=1;

% Young's Modulus
ori.panelE=2*10^9;
ori.creaseE=2*10^9;

% Initialize the size
% ori.barArea=zeros(2,1);
% ori.barLength=zeros(2,1);
% ori.barType=zeros(2,1);
% ori.sprK=zeros(2,1);


% Geometry of a single unit
ori.newNode(1,:)=[-1,0,0];
ori.newNode(2,:)=[0,0,centerSeparation/2];
ori.newNode(3,:)=[1,0,0];
ori.newNode(4,:)=[0,-1,0];
ori.newNode(5,:)=[0,0,-centerSeparation/2];
ori.newNode(6,:)=[0,1,0];

% Other initialization
newNodeNum = size(ori.newNode);
newNodeNum = newNodeNum(1);
ori.currentAppliedForce = zeros(newNodeNum,3);   
ori.currentU = zeros(newNodeNum,3);
ori.currentSprZeroStrain = pi*logical(ori.sprK);

% Add bar element
ori.AddBar(1,2,barArea,norm(ori.newNode(1,:)-ori.newNode(2,:)));
ori.AddBar(2,3,barArea,norm(ori.newNode(2,:)-ori.newNode(3,:)));
ori.AddBar(4,5,barArea,norm(ori.newNode(4,:)-ori.newNode(5,:)));
ori.AddBar(5,6,barArea,norm(ori.newNode(5,:)-ori.newNode(6,:)));

% Add Rot Spring
ori.AddHinge(1,2,5,4,barArea/100,centerSeparation,hingeStiff,pi+pi/2);
ori.AddHinge(4,2,5,3,barArea/100,centerSeparation,hingeStiff,pi+pi/2);
ori.AddHinge(3,2,5,6,barArea/100,centerSeparation,hingeStiff,pi+pi/2);
ori.AddHinge(6,2,5,1,barArea/100,centerSeparation,hingeStiff,pi+pi/2);


% Add Three Point Rot Spring
ori.threeNodeRotSpringOpen=1;
ori.threeNodeRotSpringNode(1,:)=[1,2,5];
ori.threeNodeRotSpringNode(2,:)=[3,2,5];
ori.threeNodeRotSpringNode(3,:)=[4,5,2];
ori.threeNodeRotSpringNode(4,:)=[6,5,2];

ori.threeNodeRotSpringK=hingeStiff*ones(4,1);
ori.threeNodeRotSpringTheta0=atan(5)*ones(4,1);



% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=2; % plotting range
ori.displayRangeRatio=1; % plotting range in the negative axis

ori.plotBars=1;
ori.Plot_MeshedOrigami();


%% Loading steps

% modal analysis
StiffMat=ori.Solver_CalcK();

StiffMat=(StiffMat+StiffMat')/2;
rank(StiffMat)
[mode,value]=eigs(StiffMat,7,'smallestabs');

Umode=mode(:,1);
Umode=reshape(Umode,3,6)';
ori.Plot_DeformedShape(ori.newNode,ori.newNode+Umode);

% applying the gravity loading
dc=ControllerDCLoading;

dc.increStep=100;
dc.tol=10^-4;
dc.iterMax=50;

dc.supp=[1,1,1,1;
      3,0,0,0;
      4,0,1,1;   
      6,0,1,1;];

loadMag=1;
dc.load=[3,loadMag,0,0;];
  
dc.videoOpen=1;

dc.detailFigOpen=1;
dc.selectedRefDisp=[3,1];


ori.loadingController{1}={"NR",dc};
ori.Solver_Solve();


Uselect=dc.Uhis(end,3,1);
Stiffness=loadMag*dc.increStep/Uselect;


