%% Initialize the solver
clear;clc;close all;
ori=OrigamiSolver;

%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% thickness of gold and su-8
tg=0.2*10^-6;
ts=0.80*10^-6;
% thickness of panel
tpanel=21*10^-6;
% density of SU-8
rhoSU8=1200;

% residual folding developed after fabrication (degree)
residualFold=-10;
% width of crease
W=400*10^-6;
% length of panel
Lpanel=1*10^(-3);
% undercut of XeF2 etching
underCut=100*10^-6;
% power input (mW)
qload=20;

% density of crease and panel
ori.densityCrease=rhoSU8;
ori.densityPanel=rhoSU8;

ori.node0=[0 0 0;
      Lpanel+W/2 0 0;
      2*Lpanel+W 0 0;
      0 Lpanel 0;
      Lpanel+W/2 Lpanel 0;
      2*Lpanel+W Lpanel 0;
      0*Lpanel,-0.7*Lpanel,-underCut;
      0*Lpanel,1.8*Lpanel,-underCut;
      2.5*Lpanel,-0.6*Lpanel,-underCut;
      2.5*Lpanel,1.6*Lpanel,-underCut;];
  
ori.panel0{1}=[1 2 5 4];
ori.panel0{2}=[2 3 6 5];
ori.panel0{3}=[7 8 10 9];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=3*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[500*10^(-6);tpanel;500*10^(-6)]; 
ori.panelW=W;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=1000;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(tg+ts);



%% Assign Thermal Properties

ori.panelThermalConductVec = [1.3;0.3;1.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

ori.creaseThermalCapacity = 2*10^6;
ori.panelThermalCapacity = 2*10^6;
ori.envThermalCapacity = 700;

% thickness of the submerged environment at RT
ori.t2RT=1500*10^(-6); 



%% Setup the loading controller
dynamicThermal=ControllerElectroThermalDynamicFolding;
dynamicThermal.supp=[1,1,1,1;
      2,1,1,1;
      3,1,1,1;
      4,1,1,1;
      9,1,1,1;
      10,1,1,1;
      11,1,1,1;
      12,1,1,1;];


dynamicThermal.step=80000;
dynamicThermal.dt=2*10^-6;

dynamicThermal.thermalBoundaryPanelVec=[3];
dynamicThermal.roomTempNode=[1;4];

dynamicThermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
dynamicThermal.deltaAlpha(3)=(52-14)*10^(-6); 

dynamicThermal.Emat1=39.5*10^9; 
dynamicThermal.Emat2=2*10^9;
dynamicThermal.tmat1=tg;
dynamicThermal.tmat2=ts;

% Damping
dynamicThermal.alpha=50;
dynamicThermal.beta=0.00001;

% Video Setting
dynamicThermal.videoOpen=1;
dynamicThermal.videoCropRate=200;

% the target loading of crease heating
dynamicThermal.targetCreaseHeatingHis=[3, qload/1000*ones(1,dynamicThermal.step)];

ori.loadingController{1}={"DynamicsThermal",dynamicThermal};
ori.Solver_Solve();


Uhis=dynamicThermal.Uhis;
UhisSelect=Uhis(:,6,3);
time=dynamicThermal.dt*(1:80001);
time=time';

totalTime=dynamicThermal.step*dynamicThermal.dt