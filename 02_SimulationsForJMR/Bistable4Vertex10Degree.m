%%%%%%%%%%%%%%%%%%%%%%  Origami Contact Simulator  %%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Yi Zhu, and Evegueni T. Filipov
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for nonrigid origami simulators.
% Their works pave the ground for the developement of this code.
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

a=50*10^(-3);
b=50*10^(-3);
z=10.6781*10^(-3);
Node=[0 0 z;
      -b 0 0;
      0 -b 0;
      b 0 0;
      0 b 0;
      -a -a -z;
      a -a -z;
      a a -z;
      -a a -z;];
  
Panel{1}=[2 1 5 9];
Panel{2}=[1 4 8 5];
Panel{3}=[6 3 1 2];
Panel{4}=[3 7 4 1];


u=Node(2,:)-Node(1,:);
v=Node(3,:)-Node(1,:);
CosTheta = dot(u,v)/(norm(u)*norm(v));
ThetaInDegrees = acosd(CosTheta);
MissingBeta=4*(90-ThetaInDegrees)

%% Setting up the plots for display
ViewControl=zeros(10,1);
ViewControl(1)=15; % View1: View angle 1
ViewControl(2)=30; % View2: View angle 2
ViewControl(3)=60*10^(-3); % Vsize: displayed axis range 
ViewControl(4)=1; % Vratio: ratio of displayed negative axis range versus the positive axis range

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

ratio=0.12;
RotationZeroStrain(4)=pi-ratio*pi;
TotalFoldingNum=max(FoldingSequence);
% Maximum number of loop needed for sequantial folding

%% Generate the Improved Meshing
% input parameters for generating the improved meshing
% Bar Areas, zero strain stretching will also be generated.
% Crease zero strain rotational position will be calculated.
% Crease rotational stiffness will be calculated (linear model used)

ModelConstant{1}=6*10^(-3); % CreaseW: Width of compliant creases
ModelConstant{2}=2852*10^6; % PanelE: Young's modulus of panel
ModelConstant{3}=25*10^6; % CreaseE: Young's modulus of creases
ModelConstant{4}=2.2*10^(-3); % PanelThick: thickness of panel
ModelConstant{5}=1*10^(-3); % CreaseThick: thickness of creases
ModelConstant{6}=0.4; % PanelPoisson: Poisson ratio of panel
ModelConstant{7}=0.4; % CreasePoisson: Poisson ratio of crease

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
Supp=[2,1,1,1;
      5,1,1,1;
      11,1,1,1;
      16,1,1,1;];
  
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
Supp=[2,1,1,1;
      5,1,1,1;
      11,1,1,1;
      16,1,1,1;];
  
LoadForce=10;

Load=[29,0,0,LoadForce;
      30,0,0,LoadForce;
      31,0,0,LoadForce;
      32,0,0,LoadForce;];


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
DispControler=[29,3];
ModelConstant{19}=DispControler;

[U,UhisLoading,Loadhis,StrainEnergyLoading,NodeForce,LoadForce,lockForce]...
    =NonlinearSolverLoadingDC(Panel,newNode,BarArea,BarConnect,BarLength, ...
    BarType,SprIJKL,SprK,SprTargetZeroStrain,inverseNumbering,newNumbering, ...
    Supp,Load,U,CreaseRef,CreaseNum,OldNode,LoadConstant,...
    ModelConstant,SuppElastic);
 
% [U,UhisLoading,Loadhis,StrainEnergyLoading,NodeForce,LoadForce,lockForce]...
%     =NonlinearSolverLoadingNR(Panel,newNode,BarArea,BarConnect,...
%     BarLength,BarType,SprIJKL,SprK,SprTargetZeroStrain,inverseNumbering,newNumbering, ...
%     Supp,Load,U,CreaseRef,CreaseNum,OldNode,LoadConstant,ModelConstant,SuppElastic);

%% Plotting the results
deformNode=U+newNode;
plotDeformedShape(ViewControl,AssembleNode,deformNode,newPanel,PanelNum);
% plotLoadAndReaction(ViewControl,newNode,deformNode,newPanel,Load,Supp,NodeForce,LoadForce,PanelNum);
% plotLockForce(ViewControl,newNode,deformNode,newPanel,lockForce,PanelNum);
plotLoadHis(Loadhis,UhisLoading);
plotEnergy(UhisLoading,StrainEnergyLoading,UhisAssemble,StrainEnergyAssemble);
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
