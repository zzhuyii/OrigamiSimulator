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
tic


%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% thickness of panel
tpanel=6*10^-6;

% density of polyer for the structure
rho=1200;

% length of panel (about 500um square)
Lpanel=0.5*10^(-3);
% width of crease
W=150*10^-6;

% undercut of XeF2 etching
underCut=100*10^-6;


ori.node0=[0 0 0;
      Lpanel+W/2 0 0;
      2*Lpanel+1.5*W 0 0;
      3*Lpanel+2*W 0 0;
      0 Lpanel 0;
      Lpanel+W/2 Lpanel 0;
      2*Lpanel+1.5*W Lpanel 0;
      3*Lpanel+2*W Lpanel 0;];
    
  
ori.panel0{1}=[1 2 6 5];
ori.panel0{2}=[2 3 7 6];
ori.panel0{3}=[3 4 8 7];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();

% Plot the results for inspection
ori.viewAngle1=175;
ori.viewAngle2=25;
ori.displayRange=2.5*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;

%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=W; %Length of PZT crease is 150 um
ori.creaseWidthVec(6)=W; %Length of electro-thermal crease is 250um

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
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

% Compute the equivalent thickness of folding crease
pztThick=[0.5,0.15,1,0.15]*10^-6;
pztE=[180,172,80,172]*10^9/2;
[EIpzt,C]=CalcMultiLayerBeamStiff(pztThick,pztE,1000);
TpztEquivalent=(EIpzt/ori.creaseE*12)^0.33;

etThick=[0.8,0.16]*10^-6;
etE=[2,79]*10^9/2;
[EIet,C]=CalcMultiLayerBeamStiff(etThick,etE,1000);
TetEquivalent=(EIet/ori.creaseE*12)^0.33;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=TpztEquivalent;
ori.creaseThickVec(6)=TetEquivalent;

%% setup panel contact information
ori.contactOpen=0;
ori.ke=0.0001;
ori.d0edge=40*(10^(-6));
ori.d0center=40*(10^(-6));


% %% Assign Thermal Properties
% 
% ori.panelThermalConductVec = [1.3;0.3;0.3;1.3]; 
% ori.creaseThermalConduct=0.3;
% ori.envThermalConduct=0.026;
% 
% % thickness of the submerged environment at RT
% ori.t2RT=1500*10^(-6); 


%% Frequency analysis
ori.densityCrease=rho*0.7;
ori.densityPanel=rho*0.7;

frequency=ControllerFrequencyAnalysis;
frequency.supp=[1,1,1,1;
              2,1,1,1;
              3,1,1,1;
              4,1,1,1;];
 
ori.Solver_Solve()

[frequencySquared,Umode]=ori.Dynamic_FrequencyAnalysis(frequency);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));

freq1=sqrt(freq(1))/pi/2;
freq2=sqrt(freq(2))/pi/2;

Umode1=Umode(:,index(1));
Umode1=reshape(Umode1,3,21)'/10000;

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode1);


%% Apply sine wave loading
dynamics=ControllerDynamics();
dynamics.supp=[1,1,1,1;
          2,1,1,1;
          3,1,1,1;
          4,1,1,1;];   
      


% alpha is mass; beta is stiffness
dynamics.alpha=200;
dynamics.beta=0;

% dynamics.alpha=100;
% dynamics.beta=0.00001;

omega=freq1*2*pi;
dampingRatio=dynamics.alpha/2/omega+dynamics.beta*omega/2;

% set up time steps
totalT=0.02;
dynamics.dt=1*10^-6;
step=totalT/dynamics.dt;
dynamics.Fext=zeros(step,21,3);

% Sine wave input
time=(1:step)*dynamics.dt;
dynamics.rotTargetAngle=pi*ones(step,10);
dynamics.rotTargetAngle(:,3)=1.2*pi+3/180*3.14*sin(2*pi*350*time);
% dynamics.rotTargetAngle(:,3)=1.2*pi+10/180*3.14;
dynamics.rotTargetAngle(:,6)=1.4*pi;


% ploting option
dynamics.plotOpen=0;
dynamics.videoOpen=1;
dynamics.videoCropRate=200;

% Solve the solution
selfFold=ControllerSelfFolding;
selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);

selfFold.targetRotZeroStrain(3)=1.2*pi;
selfFold.targetRotZeroStrain(6)=1.4*pi;

selfFold.supp=[1,1,1,1;
               2,1,1,1;
               3,1,1,1;
               4,1,1,1;];

selfFold.videoOpen=0;
selfFold.increStep=40;
selfFold.tol=1*10^-5;
selfFold.iterMax=50;

ori.loadingController{1}={"SelfFold",selfFold};
ori.Solver_Solve()
ori.continuingLoading=1;

ori.loadingController{1}={"Dynamics",dynamics};
ori.Solver_Solve()
toc


% tip displacement curve
Uhis=dynamics.Uhis;
node1=squeeze(Uhis(:,10,:));
node2=squeeze(Uhis(:,9,:));
node3=squeeze(Uhis(:,6,:));
node4=squeeze(Uhis(:,5,:));
v1=node2-node1;
v2=node4-node3;
angle=zeros(size(time));
for i=1:step
    v1temp=squeeze(v1(i,:));
    v2temp=squeeze(v2(i,:));
    angle(i)=acos(dot(v2temp,v1temp)/sqrt(dot(v2temp,v2temp))/sqrt(dot(v1temp,v1temp)));    
end

figure
plot(time,angle);

time=time';
angle=angle';


