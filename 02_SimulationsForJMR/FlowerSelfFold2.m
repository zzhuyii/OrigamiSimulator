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
tic
a=12*10^(-3);
b=30*10^(-3);
c=45*10^(-3);
Node=[0 0 0;
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
      
Panel{1}=[1 2 3];
Panel{2}=[1 3 4];
Panel{3}=[1 4 5];
Panel{4}=[1 5 6];
Panel{5}=[1 6 7];
Panel{6}=[1 7 8];
Panel{7}=[1 8 9];
Panel{8}=[1 9 2];

Panel{9}=[2 10 11 3];
Panel{10}=[3 12 13 4];
Panel{11}=[4 14 15 5];
Panel{12}=[5 16 17 6];
Panel{13}=[6 18 19 7];
Panel{14}=[7 20 21 8];
Panel{15}=[8 22 23 9];
Panel{16}=[9 24 25 2];

Panel{17}=[2 26 10];
Panel{18}=[3 11 27];
Panel{19}=[3 27 12];
Panel{20}=[4 13 28];
Panel{21}=[4 28 14];
Panel{22}=[5 15 29];
Panel{23}=[5 29 16];
Panel{24}=[6 17 30];
Panel{25}=[6 30 18];
Panel{26}=[7 19 31];
Panel{27}=[7 31 20];
Panel{28}=[8 21 32];
Panel{29}=[8 32 22];
Panel{30}=[9 23 33];
Panel{31}=[9 33 24];
Panel{32}=[2 25 26];


Panel{33}=[10 34 11];
Panel{34}=[12 35 13];
Panel{35}=[14 36 15];
Panel{36}=[16 37 17];
Panel{37}=[18 38 19];
Panel{38}=[20 39 21];
Panel{39}=[22 40 23];
Panel{40}=[24 41 25];

%% Setting up the plots for display
ViewControl=zeros(10,1);
ViewControl(1)=300; % View1: View angle 1
ViewControl(2)=30; % View2: View angle 2
ViewControl(3)=50*10^(-3); % Vsize: displayed axis range 
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

ratio=1;
ratio2=0.5;
RotationZeroStrain(2)=pi-ratio2*pi;
RotationZeroStrain(41)=pi-ratio2*pi;

RotationZeroStrain(8)=pi+ratio2*pi;
RotationZeroStrain(52)=pi+ratio2*pi;

RotationZeroStrain(1)=pi-ratio*pi;
RotationZeroStrain(43)=pi-ratio*pi;

RotationZeroStrain(4)=pi+ratio2*pi;
RotationZeroStrain(46)=pi+ratio2*pi;

RotationZeroStrain(6)=pi-ratio*pi;
RotationZeroStrain(49)=pi-ratio*pi;

RotationZeroStrain(10)=pi-ratio*pi;
RotationZeroStrain(55)=pi-ratio*pi;

RotationZeroStrain(12)=pi+ratio2*pi;
RotationZeroStrain(58)=pi+ratio2*pi;

RotationZeroStrain(14)=pi-ratio*pi;
RotationZeroStrain(61)=pi-ratio*pi;

TotalFoldingNum=max(FoldingSequence);
% Maximum number of loop needed for sequantial folding

%% Generate the Improved Meshing
% input parameters for generating the improved meshing
% Bar Areas, zero strain stretching will also be generated.
% Crease zero strain rotational position will be calculated.
% Crease rotational stiffness will be calculated (linear model used)

ModelConstant{1}=1.5*10^(-3); % CreaseW: Width of compliant creases
ModelConstant{2}=2000*10^6; % PanelE: Young's modulus of panel
ModelConstant{3}=20*10^6; % CreaseE: Young's modulus of creases

ModelConstant{4}=ones(40,1)*1000*10^(-6); % PanelThick: thickness of panel;

ModelConstant{5}=300*10^(-6); % CreaseThick: thickness of creases;
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
ModelConstant{13}=0.7*ModelConstant{1}; % d0center: d0 for points at the center

[newNode,newPanel,BarType,BarConnect,BarArea,BarLength,SprIJKL,...
    SprTargetZeroStrain,SprK,Type1BarNum,oldCrease,PanelInerBarStart,...
    CenterNodeStart,NewFoldingSequence,OldNode,PanelNum] ...
    =ImprovedMeshingN5B8(Node,Panel,RotationZeroStrain,...
    FoldingSequence,ModelConstant);
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
      3,1,1,1;];
Load=[1,0,0,0];
  
ModelConstant{18}=0; % NonRigidSupport
% 0 means the non rigid support is not activated.
% 1 means the non rigid support is activated.
SuppElastic=[1,3,10000;
             4,3,10000];
% first column stores node number
% second column stores direction
% third column stores stiffness

%% Nonlinear Solver for Assemble
AssembleConstant=zeros(3,1);
AssembleConstant(1)=80; % IncreStep
AssembleConstant(2)=5*10^-5; % Tor
AssembleConstant(3)=50; % iterMax

[U,UhisAssemble,StrainEnergyAssemble]=NonlinearSolverAssemble(...
    Panel,newNode,BarArea,BarConnect,BarLength,BarType,...
    SprIJKL,SprK,SprTargetZeroStrain,Supp,CreaseRef,...
    CreaseNum,NewFoldingSequence,OldNode,...
    AssembleConstant,ModelConstant,SuppElastic,Load);

AssembleNode=U+newNode;
plotDeformedShapeOnly(ViewControl,newNode,AssembleNode,newPanel,PanelNum)

%% Support and loading information for loading process
LoadForce=0.1*10^(-3);
Load=[1,0,0,-LoadForce;];

%% Nonlinear Solver for loading
LoadConstant=zeros(4,1);
LoadConstant(1)=1; % IncreStep
LoadConstant(2)=10^-8; % Tor
LoadConstant(3)=50; % iterMax
LoadConstant(4)=0.01; % LambdaBar

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
    =NonlinearSolverLoadingNR(Panel,newNode,BarArea,BarConnect,...
    BarLength,BarType,SprIJKL,SprK,SprTargetZeroStrain, ...
    Supp,Load,U,CreaseRef,CreaseNum,OldNode,LoadConstant,...
    ModelConstant,SuppElastic);

%% Plotting the results
deformNode=U+newNode;
% plotDeformedShape(ViewControl,AssembleNode,deformNode,newPanel,PanelNum);
% plotLoadAndReaction(ViewControl,newNode,deformNode,newPanel,Load,Supp,NodeForce,LoadForce,PanelNum);
plotEnergy(UhisLoading,StrainEnergyLoading,UhisAssemble,StrainEnergyAssemble);
% plotLoadHis(Loadhis,UhisLoading);
% plotLockForce(ViewControl,newNode,deformNode,newPanel,lockForce,PanelNum);

toc

AngleNum=zeros(4,1);
AngleNum(1)=543;
AngleNum(2)=544;
% Midline Crease Num
AngleNum(3)=60;
AngleNum(4)=64;
% Crease Num at the two panels
[Theta]=plotDeformedHisAndCalcTheta(ViewControl,newNode,...
    newPanel,UhisLoading,UhisAssemble,PanelNum,SprIJKL,AngleNum);


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
