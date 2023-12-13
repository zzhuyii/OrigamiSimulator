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


% Definition of Geometry
L1=300*10^-6; % Panel side overhang
L2=100*10^-6; % Beam width
L3=600*10^-6; % Panel middle width
L4=400*10^-6; % Panel length
L5=440*10^-6; % Beam Length
BeamSegment=4;
supportL=200*10^-6;


ori.node0=[L1 0 0;
           L1+L2 0 0;
           L1+L2+L3 0 0;
           L1+L2*2+L3 0 0;
           L1 supportL 0;
           L1+L2 supportL 0;
           L1+L2+L3 supportL 0;
           L1+L2*2+L3 supportL 0;];

for i=1:BeamSegment-1
    ori.node0=[ori.node0;
               L1 supportL+L5/BeamSegment*i 0;
               L1+L2 supportL+L5/BeamSegment*i 0;
               L1+L2+L3 supportL+L5/BeamSegment*i 0;
               L1+L2*2+L3 supportL+L5/BeamSegment*i 0;];
end
    ori.node0=[ori.node0;
           0 supportL+L5 0;
           L1 supportL+L5 0;
           L1+L2 supportL+L5 0;
           L1+L2+L3 supportL+L5 0;
           L1+L2*2+L3 supportL+L5 0;
           L1*2+L2*2+L3 supportL+L5 0;
           0 supportL+L4+L5 0;
           L1 supportL+L4+L5 0;
           L1+L2 supportL+L4+L5 0;
           L1+L2+L3 supportL+L4+L5 0;
           L1+L2*2+L3 supportL+L4+L5 0;
           L1*2+L2*2+L3 supportL+L4+L5 0;];

  
ori.panel0{1}=[1 2 6 5];
ori.panel0{2}=[2 3 7 6];
ori.panel0{3}=[3 4 8 7];

tempPanel=4;
for i=1:BeamSegment-1
    ori.panel0{tempPanel}=[5 6 10 9]+(i-1)*4;
    ori.panel0{tempPanel+1}=[7 8 12 11]+(i-1)*4;
    tempPanel=tempPanel+2;
end

    ori.panel0{tempPanel}=[5 6 10+1 9+1]+(BeamSegment-1)*4;
    ori.panel0{tempPanel+1}=[7 8 12+1 11+1]+(BeamSegment-1)*4;
    tempPanel=tempPanel+2;

ori.panel0{tempPanel}=[13 14 20 19]+(BeamSegment-2)*4;
ori.panel0{tempPanel+1}=[14 15 21 20]+(BeamSegment-2)*4;
ori.panel0{tempPanel+2}=[15 16 22 21]+(BeamSegment-2)*4;
ori.panel0{tempPanel+3}=[16 17 23 22]+(BeamSegment-2)*4;
ori.panel0{tempPanel+4}=[17 18 24 23]+(BeamSegment-2)*4;


% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();

% Plot the results for inspection
ori.viewAngle1=115;
ori.viewAngle2=35;
ori.displayRange=1.8*10^(-3); % plotting range
ori.displayRangeRatio=0.3; % plotting range in the negative axis

ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;

%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
creaseNumber=[3,6,4,10];
for i=1:BeamSegment
    creaseNumber=[creaseNumber,13+(i-1)*6,16+(i-1)*6];
end
creaseNumber=[creaseNumber,13+BeamSegment*6,13+BeamSegment*6+2,...
    13+BeamSegment*6+5,13+BeamSegment*6+7];


FactorForce=(BeamSegment*(BeamSegment+1)*(BeamSegment*2+1)/2/BeamSegment/BeamSegment);
FactorMoment=(BeamSegment*(BeamSegment+1)/BeamSegment);
Factor=0.5*(FactorMoment+FactorForce);

ori.creaseWidthVec(creaseNumber)=L5/Factor; 
ori.compliantCreaseOpen=0;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties
ori.panelE=150*10^9; 
ori.creaseE=70*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

% definition of panel
tpanel=2.7*10^-6;
ExtraThickFactor=14.8;

ori.panelThickVec=ExtraThickFactor*tpanel*ones(1,3); 

for i=1:BeamSegment
    ori.panelThickVec=[ori.panelThickVec,tpanel,tpanel];
end

ori.panelThickVec=[ori.panelThickVec,ExtraThickFactor*tpanel,ExtraThickFactor*tpanel,...
    ExtraThickFactor*tpanel,ExtraThickFactor*tpanel,ExtraThickFactor*tpanel];

ori.panelW=10*10^-6; 

% definition of creases
tcrease=tpanel;
ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(creaseNumber)=tcrease;
ori.creaseThickVec(1:3)=ExtraThickFactor*tcrease;
ori.creaseThickVec(end-5:end)=ExtraThickFactor*tcrease;

% density of polyer for the structure
rho=2300;


%% Frequency analysis
ori.densityCrease=rho;
ori.densityPanel=rho;

frequency=ControllerFrequencyAnalysis;
frequency.supp=[1,1,1,1;
              2,1,1,1;
              3,1,1,1;
              4,1,1,1;
              5,1,1,1;
              6,1,1,1;
              7,1,1,1;
              8,1,1,1;];
 
ori.Solver_Solve()

[frequencySquared,Umode]=ori.Dynamic_FrequencyAnalysis(frequency);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));

freq1=sqrt(freq(1))/pi/2;
freq2=sqrt(freq(2))/pi/2;
freq3=sqrt(freq(3))/pi/2;

Umode1=Umode(:,index(1));
Usize=size(Umode1,1);

Umode1=reshape(Umode1,3,Usize/3)'/10000;

Umode2=Umode(:,index(2));
Umode2=reshape(Umode2,3,Usize/3)'/10000;

Umode3=Umode(:,index(3));
Umode3=reshape(Umode3,3,Usize/3)'/10000;

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode1);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode2);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode3);



%% Statically Fold to an Angle
selfFold=ControllerSelfFolding();
foldAngleRate=1.05;

selfFold.supp=[1,1,1,1;
              2,1,1,1;
              3,1,1,1;
              4,1,1,1;
              5,1,1,1;
              6,1,1,1;
              7,1,1,1;
              8,1,1,1;]; 

activeCreaseNumberAll=[4,10];
for i=1:BeamSegment
    activeCreaseNumberAll=[activeCreaseNumberAll,13+(i-1)*6,16+(i-1)*6];
end

selfFold.increStep=20;
selfFold.tol=1*10^-5;

selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);
selfFold.targetRotZeroStrain(activeCreaseNumberAll)=foldAngleRate*pi;


selfFold.videoOpen=0;
selfFold.plotOpen=1;

ori.loadingController{1}={"SelfFold",selfFold};
ori.Solver_Solve();

ori.continuingLoading=1;

%% Solve Frequency

[frequencySquared,Umode]=ori.Dynamic_FrequencyAnalysis(frequency);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));

freq1Fold=sqrt(freq(1))/pi/2;
freq2Fold=sqrt(freq(2))/pi/2;
freq3Fold=sqrt(freq(3))/pi/2;

Umode1Fold=Umode(:,index(1));
Usize=size(Umode1Fold,1);

Umode1Fold=reshape(Umode1Fold,3,Usize/3)'/10000;

Umode2Fold=Umode(:,index(2));
Umode2Fold=reshape(Umode2Fold,3,Usize/3)'/10000;

Umode3Fold=Umode(:,index(3));
Umode3Fold=reshape(Umode3Fold,3,Usize/3)'/10000;

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode1Fold);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode2Fold);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode3Fold);

      

%% Apply sine wave loading
dynamics=ControllerDynamics();
dynamics.supp=[1,1,1,1;
              2,1,1,1;
              3,1,1,1;
              4,1,1,1;
              5,1,1,1;
              6,1,1,1;
              7,1,1,1;
              8,1,1,1;]; 
      

% alpha is mass; beta is stiffness
% this is for Rayleigh Damping

dynamics.alpha=200;
dynamics.beta=0;
omega=freq1*2*pi;
dampingRatio=dynamics.alpha/2/omega+dynamics.beta*omega/2;


% Set up loading time
dynamics.dt=5*10^-7;
totalTime=0.02;
step=round(totalTime/dynamics.dt);

% External forces
dynamics.Fext=zeros(step,Usize/3,3); 

% Sine wave input
time=(1:step)*dynamics.dt;
dynamics.rotTargetAngle=pi*ones(step,ori.oldCreaseNum);


Freq=freq1;
activeCreaseNumber=[4];
for i=1:BeamSegment
    activeCreaseNumber=[activeCreaseNumber,13+(i-1)*6];
end

activeCreaseNumber2=[10];
for i=1:BeamSegment
    activeCreaseNumber2=[activeCreaseNumber2,16+(i-1)*6];
end

for i=1:length(activeCreaseNumber)
    dynamics.rotTargetAngle(:,activeCreaseNumber(i))=(foldAngleRate*pi+1/180*3.14*sin(2*pi*Freq*time))';
end

for i=1:length(activeCreaseNumber2)
    dynamics.rotTargetAngle(:,activeCreaseNumber2(i))=(foldAngleRate*pi+1/180*3.14*sin(2*pi*Freq*time))';
end


% ploting option
dynamics.plotOpen=1;
dynamics.videoOpen=0;
dynamics.videoCropRate=100;


ori.loadingController{1}={"Dynamics",dynamics};
ori.Solver_Solve()
toc
 
% ori.Plot_DeformedHis(ori.newNode,dynamics.Uhis(20000:100:end,:,:))


% panel rotation
Uhis=dynamics.Uhis;

angle=zeros(step,1);
for i=1:step
    node1=squeeze(Uhis(i,32,:))+ori.newNode(32,:)';
    node2=squeeze(Uhis(i,26,:))+ori.newNode(26,:)';
    vector=node1-node2;    
    angle(i)=atan(vector(3)/vector(2))*180/pi;
    if vector(2)<0
        angle(i)=atan(vector(3)/vector(2))*180/pi+180;
    end
end

figure
hold on
plot(time,angle)

UhisCorp=Uhis(1:100:end,:,:);
ori.Plot_DeformedHis(ori.newNode,UhisCorp);


