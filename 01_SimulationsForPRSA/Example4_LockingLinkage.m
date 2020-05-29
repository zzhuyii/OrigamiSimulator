%%%%%%%%%%%%%%%%%%%%%%  Origami Contact Simulator  %%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Yi Zhu, and Evegueni T. Filipov
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% contact model presented in this simulation code. 
%
% Code Features: 
% (1) Simulate Contact in Origami patterns
% (2) Simulate Compliant Creases in Origami Patterns
% (3) Automated Meshing for Compliant Creases
% (4) Solver for both folding and loading available
% (5) Nonrigid support available
%
% Reference:
% [1] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' PRSA. (submitted)       
% [2] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [3] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  
% [4] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [5] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%%%%%%%%%%%%%%%%%%%  Origami Contact Simulator  %%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;
%% Input data for Defining the Geometry
% This section of code is used to generate the geometry of structure
% Input files are still origami pattern with concentrated hinges.The code
% will automatically generate compliant creases.

a=20*10^(-3);
b=4*10^(-3);
c=10*10^(-3);
z=3*10^(-3);
Node=[b 0 0;
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
  
Panel{1}=[1 2 8 7];
Panel{2}=[2 3 9 8];
Panel{3}=[3 4 10 9];
Panel{4}=[4 5 11 10];
Panel{5}=[5 6 12 11];
Panel{6}=[7 8 14 13];
Panel{7}=[8 9 15 14];
Panel{8}=[9 10 16 15];
Panel{9}=[10 11 17 16];
Panel{10}=[11 12 18 17];

Panel{11}=[1 2 8 7]+[18 18 18 18];
Panel{12}=[2 3 9 8]+[18 18 18 18];
Panel{13}=[3 4 10 9]+[18 18 18 18];
Panel{14}=[4 5 11 10]+[18 18 18 18];
Panel{15}=[5 6 12 11]+[18 18 18 18];
Panel{16}=[7 8 14 13]+[18 18 18 18];
Panel{17}=[8 9 15 14]+[18 18 18 18];
Panel{18}=[9 10 16 15]+[18 18 18 18];
Panel{19}=[10 11 17 16]+[18 18 18 18];
Panel{20}=[11 12 18 17]+[18 18 18 18];


%% Setting up the plots for display
ViewControl=zeros(10,1);
ViewControl(1)=15; % View1: View angle 1
ViewControl(2)=15; % View2: View angle 2
ViewControl(3)=70*10^(-3); % Vsize: displayed axis range 
ViewControl(4)=0.1; % Vratio: ratio of displayed negative axis range versus the positive axis range

%% Assign Zero strain position for creases
[CreaseNum,Crease,CreaseType]=IdentifyCrease(Node,Panel);
% Here we identify the creases and show the drawing of the creases. With
% this information, users can assign mountain and valley folds and their
% zero strain position manually if needed.
plotOriginalMeshing(Node,Panel,CreaseNum,Crease,ViewControl)

RotationZeroStrain=pi*ones(CreaseNum,1);
% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
FoldingSequence=ones(CreaseNum,1);
% Folding Sequence indicate which crease will be folded first

ratio=0.5;
RotationZeroStrain(4)=pi+ratio*pi;
RotationZeroStrain(7)=pi-ratio*pi;
RotationZeroStrain(10)=pi+ratio*pi;
RotationZeroStrain(13)=pi-ratio*pi;
RotationZeroStrain(16)=pi+ratio*pi;

RotationZeroStrain(3)=pi-ratio*pi;
RotationZeroStrain(18)=pi-ratio*pi;
RotationZeroStrain(6)=pi+ratio*pi;
RotationZeroStrain(20)=pi+ratio*pi;
RotationZeroStrain(9)=pi+ratio*pi;
RotationZeroStrain(22)=pi+ratio*pi;
RotationZeroStrain(12)=pi-ratio*pi;
RotationZeroStrain(24)=pi-ratio*pi;

RotationZeroStrain(4+27)=pi+ratio*pi;
RotationZeroStrain(7+27)=pi-ratio*pi;
RotationZeroStrain(10+27)=pi+ratio*pi;
RotationZeroStrain(13+27)=pi-ratio*pi;
RotationZeroStrain(16+27)=pi+ratio*pi;

RotationZeroStrain(3+27)=pi-ratio*pi;
RotationZeroStrain(18+27)=pi-ratio*pi;
RotationZeroStrain(6+27)=pi+ratio*pi;
RotationZeroStrain(20+27)=pi+ratio*pi;
RotationZeroStrain(9+27)=pi+ratio*pi;
RotationZeroStrain(22+27)=pi+ratio*pi;
RotationZeroStrain(12+27)=pi-ratio*pi;
RotationZeroStrain(24+27)=pi-ratio*pi;

TotalFoldingNum=max(FoldingSequence);
% Maximum number of loop needed for sequantial folding

%% Generate the Improved Meshing
% input parameters for generating the improved meshing
% Bar Areas, zero strain stretching will also be generated.
% Crease zero strain rotational position will be calculated.
% Crease rotational stiffness will be calculated (linear model used)

ModelConstant{1}=2*10^(-3); % CreaseW: Width of compliant creases
ModelConstant{2}=2*10^9; % PanelE: Young's modulus of panel
ModelConstant{3}=2*10^9; % CreaseE: Young's modulus of creases
ModelConstant{4}=500*10^(-6); % PanelThick: thickness of panel;
ModelConstant{5}=100*10^(-6); % CreaseThick: thickness of creases;
ModelConstant{6}=0.3; % PanelPoisson: Poisson ratio of panel
ModelConstant{7}=0.3; % CreasePoisson: Poisson ratio of crease

ModelConstant{8}=3; % Flag2D3D: 
% Flag2D3D is used to determine how the additional crease structuer is
% generated,2D means center bars are genertaed through using an averaged 
% vector, 3D means the center bars are located at the original positoin.
% 3 3D, 2 2D

ModelConstant{9}=4; % DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.

ModelConstant{17}=0; % CompliantCreaseOpen
% 1 means include compliant crease model
% 0 means using concentrated hinge model

ModelConstant{10}=0; % LockingOpen:
% 1: calculating the locking forces and formulate stiffness matrix
% 0: Close the calculation for locking induced by having panel interaction

ModelConstant{11}=0.02; % ke: used to scale the magnitude of potentil
ModelConstant{12}=1*(10^(-3)); % d0edge: d0 for points at the edge
ModelConstant{13}=1*(10^(-3)); % d0center: d0 for points at the center

[newNode,newPanel,BarType,BarConnect,BarArea,BarLength,SprIJKL,SprTargetZeroStrain, ... 
    SprK,Type1BarNum,oldCrease,PanelInerBarStart,CenterNodeStart,NewFoldingSequence,OldNode,PanelNum] ...
    =ImprovedMeshingN5B8(Node,Panel,RotationZeroStrain,FoldingSequence,ModelConstant);
% Generate Improved meshing

plotImprovedMeshing(ViewControl,newNode,newPanel,BarArea,BarConnect);
% Plot improved meshing for inspection
[newNumbering,inverseNumbering]=NewSequence(newNode);
% Generate newNumbering and inverseNumbering code for sparse matrix
[CreaseRef]= NumberingForCreaseLocking(oldCrease,CreaseNum,BarType);

ModelConstant{14}=TotalFoldingNum; % TotalFoldingNum
ModelConstant{15}=PanelInerBarStart; % TotalFoldingNum
ModelConstant{16}=CenterNodeStart; % TotalFoldingNum

%% Input information of support for Assemble
Supp=[1,1,1,1;
      2,1,1,1;
      13,0,0,1;
      14,0,0,1;
      19,1,1,1;
      20,1,1,1;
      31,0,0,1;
      32,0,0,1;];
  
ModelConstant{18}=0; % NonRigidSupport
ksup=500;
SuppElastic=[19,3,ksup;
             20,3,ksup;
             31,3,ksup;
             32,3,ksup;
             23,3,ksup;
             24,3,ksup;
             35,3,ksup;
             36,3,ksup;];
% first column stores node number
% second column stores direction
% third column stores stiffness

%% Nonlinear Solver for Assemble
AssembleConstant=zeros(3,1);
AssembleConstant(1)=30; % IncreStep
AssembleConstant(2)=10^-6; % Tor
AssembleConstant(3)=50; % iterMax

[U,UhisAssemble,StrainEnergyAssemble]=NonlinearSolverAssemble(...
    Panel,newNode,BarArea,BarConnect,BarLength,BarType,SprIJKL,SprK,...
    SprTargetZeroStrain,inverseNumbering,newNumbering, ...
    Supp,CreaseRef,CreaseNum,NewFoldingSequence,OldNode,...
    AssembleConstant,ModelConstant,SuppElastic);

AssembleNode=U+newNode;
plotDeformedShapeOnly(ViewControl,newNode,AssembleNode,newPanel,PanelNum)

%% Support and loading information for loading process

Supp=[2,1,1,1;
      14,1,1,1;
      5,1,1,1;
      17,1,1,1;
      32,1,1,0;
      20,1,1,0;
      23,1,1,0;
      35,1,1,0;];
  
ModelConstant{18}=1; % NonRigidSupport
ksup=5;
SuppElastic=[19,3,ksup;
             20,3,ksup;
             31,3,ksup;
             32,3,ksup;
             23,3,ksup;
             24,3,ksup;
             35,3,ksup;
             36,3,ksup;];
% first column stores node number
% second column stores direction
% third column stores stiffness

LoadForce=5*10^(-3);
Load=[21,0,0,LoadForce;
      22,0,0,LoadForce;
      33,0,0,LoadForce;
      34,0,0,LoadForce;]; 
  
ModelConstant{10}=1; % LockingOpen

%% Nonlinear Solver for loading
LoadConstant=zeros(4,1);
LoadConstant(1)=30; % IncreStep
LoadConstant(2)=10^-6; % Tor
LoadConstant(3)=50; % iterMax
LoadConstant(4)=1; % LambdaBar

% Use Theta0 to set the zero strain position of crease so that creases are
% at zero strain state at the start of loading. Activate the following code
% if no self assemble process is used.

% AssembleNode=newNode;
% A=size(newNode);
% UhisAssemble=zeros(1,A(1),A(2));
% U=zeros(size(newNode));
% StrainEnergyAssemble=zeros(1,4);
% [Theta0]=CreaseTheta(U,CreaseIJKL,newNode);

[U,UhisLoading,Loadhis,StrainEnergyLoading,NodeForce,LoadForce,lockForce]...
    =NonlinearSolverLoadingNR(Panel,newNode,BarArea,BarConnect,BarLength, ...
    BarType,SprIJKL,SprK,SprTargetZeroStrain,inverseNumbering,newNumbering, ...
    Supp,Load,U,CreaseRef,CreaseNum,OldNode,LoadConstant,...
    ModelConstant,SuppElastic);


%% Plotting the results
deformNode=U+newNode;
%plotDeformedShapeOnly(ViewControl,deformNode,deformNode,newPanel,PanelNum);
%plotLoadAndReaction(ViewControl,newNode,deformNode,newPanel,Load,Supp,NodeForce,LoadForce,PanelNum);
% plotLockForce(ViewControl,newNode,deformNode,newPanel,lockForce,PanelNum);
% plotEnergy(UhisLoading,StrainEnergyLoading,UhisAssemble,StrainEnergyAssemble);
plotDeformedHis(ViewControl,newNode,newPanel,UhisLoading,UhisAssemble,PanelNum);
% plotLoadHis(Loadhis,UhisLoading);


%% Summary of ModelConstant
% CreaseW=ModelConstant{1};
% PanelE=ModelConstant{2};
% CreaseE=ModelConstant{3};
% PanelThick=ModelConstant{4};
% CreaseThick=ModelConstant{5};
% PanelPoisson=ModelConstant{6};
% CreasePoisson=ModelConstant{7};
% Flag2D3D=ModelConstant{8};
% DiagonalRate=ModelConstant{9};
% LockingOpen=ModelConstant{10};
% ke=ModelConstant{11};
% d0edge=ModelConstant{12};
% d0center=ModelConstant{13};
% TotalFoldingNum=ModelConstant{14};
% PanelInerBarStart=ModelConstant{15};
% CenterNodeStart=ModelConstant{16};
% CompliantCreaseOpen=ModelConstant{17};
% ElasticSupportOpen=ModelConstant{18};
