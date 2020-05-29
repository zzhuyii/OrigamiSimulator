%% Nonlinear Solver for Improved Meshing
clear;clc;close all;
%% Input data for Defining the Geometry
% This section of code is used to generate the geometry of structure
% The meshing generated in this section is still based on simple meshing
% and is still using the infinitesimal hinge assumption.
tic
a=20*10^(-3);
b=2*10^(-3);
Node=[0 0 0;
      0 a 0;
      0 2*a 0;
      0 3*a 0;
      0 4*a 0;
      0 5*a 0;
      a 0 0;
      a a 0;
      a 4*a 0;
      a 5*a 0;
      2*a 0 0;
      2*a 5*a 0;
      3*a 0 0;
      3*a 5*a 0;
      4*a 0 0;
      4*a a 0;
      4*a 4*a 0;
      4*a 5*a 0;
      5*a 0 0;
      5*a a 0;
      5*a 2*a 0;
      5*a 3*a 0;
      5*a 4*a 0;
      5*a 5*a 0;];
  
 Node(:,3)=Node(:,3)+2*a*ones(24,1);
  
Panel{1}=[1 8 2];
Panel{2}=[2 8 3];
Panel{3}=[8 9 4 3];
Panel{4}=[4 9 5];
Panel{5}=[5 9 6];
Panel{6}=[7 8 1];
Panel{7}=[9 10 6];
Panel{8}=[7 11 8];
Panel{9}=[8 16 17 9];
Panel{10}=[9 12 10];
Panel{11}=[11 13 16 8];
Panel{12}=[9 17 14 12];
Panel{13}=[13 15 16];
Panel{14}=[14 17 18];
Panel{15}=[15 19 16];
Panel{16}=[19 20 16];
Panel{17}=[20 21 16];
Panel{18}=[16 21 22 17];
Panel{19}=[22 23 17];
Panel{20}=[23 24 17];
Panel{21}=[17 24 18];

% this will generate a Miura sheet model

%% Setting up the plots for display
ViewControl=zeros(10,1);
ViewControl(1)=30; % View1: View angle 1
ViewControl(2)=15; % View2: View angle 2
ViewControl(3)=100*10^(-3); % Vsize: displayed axis range 
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
ratio=0;
RotationZeroStrain(19)=pi-ratio*pi;
RotationZeroStrain(6)=pi-ratio*pi;
RotationZeroStrain(21)=pi-ratio*pi;
RotationZeroStrain(20)=pi-ratio*pi;

RotationZeroStrain(5)=pi+ratio*pi;
RotationZeroStrain(18)=pi+ratio*pi;
RotationZeroStrain(7)=pi+ratio*pi;
RotationZeroStrain(22)=pi+ratio*pi;
RotationZeroStrain(26)=pi+ratio*pi;
RotationZeroStrain(39)=pi+ratio*pi;
RotationZeroStrain(25)=pi+ratio*pi;
RotationZeroStrain(37)=pi+ratio*pi;

RotationZeroStrain(2)=pi-ratio*pi;
RotationZeroStrain(3)=pi-ratio*pi;
RotationZeroStrain(14)=pi-ratio*pi;

RotationZeroStrain(29)=pi-ratio*pi;
RotationZeroStrain(33)=pi-ratio*pi;
RotationZeroStrain(35)=pi-ratio*pi;

RotationZeroStrain(31)=pi-ratio*pi;
RotationZeroStrain(41)=pi-ratio*pi;
RotationZeroStrain(43)=pi-ratio*pi;

RotationZeroStrain(10)=pi-ratio*pi;
RotationZeroStrain(12)=pi-ratio*pi;
RotationZeroStrain(15)=pi-ratio*pi;


FoldingSequence=ones(CreaseNum,1);
% Folding Sequence indicate which crease will be folded first

ratio=0.12;
RotationZeroStrain(4)=pi-ratio*pi;
TotalFoldingNum=max(FoldingSequence);
% Maximum number of loop needed for sequantial folding

%% Generate the Improved Meshing
% input parameters for generating the improved meshing
% Bar Areas, zero strain stretching will also be generated.
% Crease zero strain rotational position will be calculated.
% Crease rotational stiffness will be calculated (linear model used)

ModelConstant{1}=4*10^(-3); % CreaseW: Width of compliant creases
ModelConstant{2}=2*10^9; % PanelE: Young's modulus of panel
ModelConstant{3}=2*10^9; % CreaseE: Young's modulus of creases
ModelConstant{4}=500*10^(-6); % PanelThick: thickness of panel;
ModelConstant{5}=100*10^(-6); % CreaseThick: thickness of creases;
ModelConstant{6}=0.3; % PanelPoisson: Poisson ratio of panel
ModelConstant{7}=0.3; % CreasePoisson: Poisson ratio of crease

ModelConstant{8}=2; % Flag2D3D: 
% Flag2D3D is used to determine how the additional crease structuer is
% generated,2D means center bars are genertaed through using an averaged 
% vector, 3D means the center bars are located at the original positoin.
% 3 3D, 2 2D

ModelConstant{9}=4; % DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.

ModelConstant{17}=1; % CompliantCreaseOpen
% 1 means include compliant crease model
% 0 means using concentrated hinge model

ModelConstant{10}= 0; % LockingOpen:
% 1: calculating the locking forces and formulate stiffness matrix
% 0: Close the calculation for locking induced by having panel interaction

ModelConstant{11}=0.002; % ke: used to scale the magnitude of potentil
ModelConstant{12}=ModelConstant{1}; % d0edge: d0 for points at the edge
ModelConstant{13}=ModelConstant{1}; % d0center: d0 for points at the center

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
Supp=[26,1,1,1;
      27,1,1,1;
      28,1,1,1;
      29,1,1,1;];
  
ModelConstant{18}=0; % NonRigidSupport
% 0 means the non rigid support is not activated.
% 1 means the non rigid support is activated.

SuppElastic=[1,3,10000;
             4,3,10000];
% first column stores node number
% second column stores direction
% third column stores stiffness

%% Nonlinear Solver for Assemble
% AssembleConstant=zeros(3,1);
% AssembleConstant(1)=20; % IncreStep
% AssembleConstant(2)=10^-6; % Tor
% AssembleConstant(3)=50; % iterMax
% 
% [U,UhisAssemble,StrainEnergyAssemble]=NonlinearSolverAssemble(...
%     Panel,newNode,BarArea,BarConnect,BarLength,BarType,SprIJKL,SprK,...
%     SprTargetZeroStrain,inverseNumbering,newNumbering, ...
%     Supp,CreaseRef,CreaseNum,NewFoldingSequence,OldNode,...
%     AssembleConstant,ModelConstant,SuppElastic);
% 
% AssembleNode=U+newNode;

% plotDeformedShapeOnly(ViewControl,newNode,AssembleNode,newPanel,PanelNum)


%% Support and loading information for loading process
Supp=[26,1,1,1;
      27,1,1,1;
      28,1,1,1;
      29,1,1,1;];
  

LoadForce=2.0*10^(-3);

Load=[33,0,LoadForce,-LoadForce;
      34,0,LoadForce,-LoadForce;
      57,-LoadForce,0,-LoadForce;
      58,-LoadForce,0,-LoadForce;
      39,0,-LoadForce,-LoadForce;
      40,0,-LoadForce,-LoadForce;
      09,LoadForce,0,-LoadForce;
      10,LoadForce,0,-LoadForce;];


%% Nonlinear Solver for loading
LoadConstant=zeros(4,1);
LoadConstant(1)=40; % IncreStep
LoadConstant(2)=10^-6; % Tor
LoadConstant(3)=50; % iterMax
LoadConstant(4)=1; % LambdaBar

% Use Theta0 to set the zero strain position of crease so that creases are
% at zero strain state at the start of loading. Activate the following code
% if no self assemble process is used.

AssembleNode=newNode;
A=size(newNode);
UhisAssemble=zeros(1,A(1),A(2));
U=zeros(size(newNode));
StrainEnergyAssemble=zeros(1,4);
[SprTargetZeroStrain]=CreaseTheta(U,SprIJKL,newNode);

% if displacement controled method is used, use the following to set up the
% controlling displacement entry
DispControler=[22,3];
ModelConstant{19}=DispControler;

[U,UhisLoading,Loadhis,StrainEnergyLoading,NodeForce,LoadForce,lockForce]...
    =NonlinearSolverLoadingMGDCM(Panel,newNode,BarArea,BarConnect,BarLength, ...
    BarType,SprIJKL,SprK,SprTargetZeroStrain,inverseNumbering,newNumbering, ...
    Supp,Load,U,CreaseRef,CreaseNum,OldNode,LoadConstant,...
    ModelConstant,SuppElastic);

%% Plotting the results
deformNode=U+newNode;
plotDeformedShape(ViewControl,AssembleNode,deformNode,newPanel,PanelNum);
% plotLoadAndReaction(ViewControl,newNode,deformNode,newPanel,Load,Supp,NodeForce,LoadForce,PanelNum);
% plotLockForce(ViewControl,newNode,deformNode,newPanel,lockForce,PanelNum);
plotLoadHis(Loadhis,UhisLoading);
plotEnergy(UhisLoading,StrainEnergyLoading,UhisAssemble,StrainEnergyAssemble);
toc
plotDeformedHis(ViewControl,newNode,newPanel,UhisLoading,UhisAssemble,PanelNum);


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


