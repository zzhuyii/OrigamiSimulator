%% Nonlinear Solver for Improved Meshing
clear;clc;close all;
%% Input data for Defining the Geometry
% This section of code is used to generate the geometry of structure
% The meshing generated in this section is still based on simple meshing
% and is still using the infinitesimal hinge assumption.
tic
a=20*10^(-3);
b=2*10^(-3);
node0=[0 0 0;
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
  
node0(:,3)=node0(:,3)+2*a*ones(24,1);
  
panel0{1}=[1 8 2];
panel0{2}=[2 8 3];
panel0{3}=[8 9 4 3];
panel0{4}=[4 9 5];
panel0{5}=[5 9 6];
panel0{6}=[7 8 1];
panel0{7}=[9 10 6];
panel0{8}=[7 11 8];
panel0{9}=[8 16 17 9];
panel0{10}=[9 12 10];
panel0{11}=[11 13 16 8];
panel0{12}=[9 17 14 12];
panel0{13}=[13 15 16];
panel0{14}=[14 17 18];
panel0{15}=[15 19 16];
panel0{16}=[19 20 16];
panel0{17}=[20 21 16];
panel0{18}=[16 21 22 17];
panel0{19}=[22 23 17];
panel0{20}=[23 24 17];
panel0{21}=[17 24 18];

% this will generate a Miura sheet model

%% Setting up the plots for display

viewControl=zeros(10,1);

% View1: View angle 1
viewControl(1)=30; 

% View2: View angle 2
viewControl(2)=15; 

% Vsize: displayed axis range 
viewControl(3)=100*10^(-3); 

% Vratio: ratio of displayed negative axis range versus the positive axis range
viewControl(4)=0.1; 


%% Generate the geometry of compliant crease

% Here we identify the creases based on the original input;
[oldCreaseNum,oldCreaseConnect,oldCreaseType]=Mesh_IdentifyCrease(node0,panel0);

% Plot the original meshing for inspection;
Plot_OriginalMeshing(node0,panel0,oldCreaseNum,oldCreaseConnect,viewControl)

% flag2D3D is used to determine how the additional crease structuer is
% generated,2D means center bars are genertaed through using an averaged 
% vector, 3D means the center bars are located at the original positoin.
% 3 3D, 2 2D
modelGeometryConstant{1}=2; 

% 1 means include compliant crease model
% 0 means using concentrated hinge model
modelGeometryConstant{2}=1; 

% crease width to generate the topology of the compliant 
% crease origami;
creaseWidthMat=zeros(oldCreaseNum,1);

creaseWidthMat(19)=4*10^(-3);
creaseWidthMat(6)=4*10^(-3);
creaseWidthMat(21)=4*10^(-3);
creaseWidthMat(20)=4*10^(-3);

creaseWidthMat(5)=4*10^(-3);
creaseWidthMat(18)=4*10^(-3);
creaseWidthMat(7)=4*10^(-3);
creaseWidthMat(22)=4*10^(-3);
creaseWidthMat(26)=4*10^(-3);
creaseWidthMat(39)=4*10^(-3);
creaseWidthMat(25)=4*10^(-3);
creaseWidthMat(37)=4*10^(-3);

creaseWidthMat(2)=4*10^(-3);
creaseWidthMat(3)=4*10^(-3);
creaseWidthMat(14)=4*10^(-3);

creaseWidthMat(29)=4*10^(-3);
creaseWidthMat(33)=4*10^(-3);
creaseWidthMat(35)=4*10^(-3);

creaseWidthMat(31)=4*10^(-3);
creaseWidthMat(41)=4*10^(-3);
creaseWidthMat(43)=4*10^(-3);

creaseWidthMat(10)=4*10^(-3);
creaseWidthMat(12)=4*10^(-3);
creaseWidthMat(15)=4*10^(-3);

modelGeometryConstant{3}=creaseWidthMat; 

% generate the geometry of system
[newNode,newPanel,barType,barConnect,...
    sprIJKL,type1BarNum,panelInnerBarStart,centerNodeStart,...
    newNode2OldNode,newCrease2OldCrease,newPanel2OldPanel,newPanelNum] ...
    =Mesh_CompliantCreaseGeometry(node0,panel0,...
    oldCreaseNum,oldCreaseConnect,oldCreaseType,...
    modelGeometryConstant);

% Plot the pattern with updated geometry
Plot_ImprovedMeshing(viewControl,newNode,newPanel,barConnect);

% generate creaseRef matrix used for calculating the contact
[creaseRef]= Mesh_NumberingForContact(newCrease2OldCrease,oldCreaseNum);

% calcualte the barLength 
barLength=Mesh_BarLength(newNode,barConnect);

% Update topology constant after generating the new topology
modelGeometryConstant{4}=panelInnerBarStart; 
modelGeometryConstant{5}=centerNodeStart; 
modelGeometryConstant{6}=type1BarNum;


%% Assign Mechanical Properties

% PanelE: Young's modulus of panel
modelMechanicalConstant{1}=2000*10^6; 

% CreaseE: Young's modulus of creases
modelMechanicalConstant{2}=2000*10^6; 

% PanelPoisson: Poisson ratio of panel
modelMechanicalConstant{3}=0.3; 

% CreasePoisson: Poisson ratio of crease
modelMechanicalConstant{4}=0.3; 

% PanelThick: thickness of panel;
% This is a vector storeing the thicknes of panels. 
modelMechanicalConstant{5}=ones(21,1)*500*10^(-6);

% thickness of creases;
creaesThicknessMat=zeros(oldCreaseNum,1);

creaesThicknessMat(19)=100*10^(-6);
creaesThicknessMat(6)=100*10^(-6);
creaesThicknessMat(21)=100*10^(-6);
creaesThicknessMat(20)=100*10^(-6);

creaesThicknessMat(5)=100*10^(-6);
creaesThicknessMat(18)=100*10^(-6);
creaesThicknessMat(7)=100*10^(-6);
creaesThicknessMat(22)=100*10^(-6);
creaesThicknessMat(26)=100*10^(-6);
creaesThicknessMat(39)=100*10^(-6);
creaesThicknessMat(25)=100*10^(-6);
creaesThicknessMat(37)=100*10^(-6);

creaesThicknessMat(2)=100*10^(-6);
creaesThicknessMat(3)=100*10^(-6);
creaesThicknessMat(14)=100*10^(-6);

creaesThicknessMat(29)=100*10^(-6);
creaesThicknessMat(33)=100*10^(-6);
creaesThicknessMat(35)=100*10^(-6);

creaesThicknessMat(31)=100*10^(-6);
creaesThicknessMat(41)=100*10^(-6);
creaesThicknessMat(43)=100*10^(-6);

creaesThicknessMat(10)=100*10^(-6);
creaesThicknessMat(12)=100*10^(-6);
creaesThicknessMat(15)=100*10^(-6);

modelMechanicalConstant{6}=creaesThicknessMat; 


% DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.
modelMechanicalConstant{7}=4; 

% panelW
% this is an averaged creaseW. used to calculate panel bending stiffness
panelW=4*10^-3;
modelMechanicalConstant{8}=panelW;



%% panel contact related input

% contact open
% 1: means consider panel contact;
% 0: means ignore panel contact
modelMechanicalConstant{9}= 0; 

% ke: used to scale the magnitude of potentil
modelMechanicalConstant{10}=0.002; 

% d0edge: d0 for points at the edge
modelMechanicalConstant{11}=6*10^-3;

% d0center: d0 for points at the center
modelMechanicalConstant{12}=6*10^-3;



%% Assign zero strain position for creases during self-folding

% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
rotationZeroStrain=pi*ones(oldCreaseNum,1);

% Folding Sequence indicate which crease will be folded first
foldingSequence=ones(oldCreaseNum,1);


% Maximum number of loop needed for sequantial folding
totalFoldingNum=max(foldingSequence);


% set the self-folding related proeprties
modelMechanicalConstant{13}=rotationZeroStrain;
modelMechanicalConstant{14}=totalFoldingNum; 
modelMechanicalConstant{15}=foldingSequence;

% distribution Factor of rotation zero strain
modelMechanicalConstant{16}=0.5; 


%% Generate mechanical properties of the origami
[barArea,sprK,sprTargetZeroStrain,sprFoldingSequence]...
    =Mesh_MechanicalProperty(modelMechanicalConstant,...
    modelGeometryConstant,oldCreaseType,...
    oldCreaseNum,creaseRef,barLength,panel0,...
    barConnect,newNode);



%% Input information of support

% define support information
supp=[26,1,1,1;
      27,1,1,1;
      28,1,1,1;
      29,1,1,1;];
supportInfo{1}=supp; 
  
% NonRigidSupport
% 0 means the non rigid support is not activated.
% 1 means the non rigid support is activated.
supportInfo{2}=0; 

% first column stores node number
% second column stores direction
% third column stores stiffness
suppElastic=[1,3,10000;
             4,3,10000];
supportInfo{3}=suppElastic;



%% Support and loading information for loading process

loadForce=2.0*10^(-3);

load=[33,0,loadForce,-loadForce;
      34,0,loadForce,-loadForce;
      57,-loadForce,0,-loadForce;
      58,-loadForce,0,-loadForce;
      39,0,-loadForce,-loadForce;
      40,0,-loadForce,-loadForce;
      09,loadForce,0,-loadForce;
      10,loadForce,0,-loadForce;];


%% Nonlinear Solver for loading

loadConstant=zeros(4,1);

% increStep
loadConstant(1)=40; 

% tor
loadConstant(2)=10^-6; 

% iterMax
loadConstant(3)=50;

% lambdaBar
loadConstant(4)=1; 


assembleNode=newNode;
A=size(newNode);
UhisAssemble=zeros(1,A(1),A(2));
U=zeros(size(newNode));
strainEnergyAssemble=zeros(1,4);
[sprTargetZeroStrain]=Spr_Theta(U,sprIJKL,newNode);


[U,UhisLoading,loadHis,strainEnergyLoading,nodeForce,loadForce,contactForce]...
    =Solver_LoadingMGDCM(panel0,newNode2OldNode,oldCreaseNum,...
    newNode,barConnect,barType,barArea,barLength,...
    sprIJKL,sprK,sprTargetZeroStrain,creaseRef,load,supportInfo,U,...
    loadConstant,modelGeometryConstant,modelMechanicalConstant);


%% Plotting the results
deformNode=U+newNode;
Plot_DeformedShape(viewControl,assembleNode,deformNode,newPanel);
Plot_LoadHis(loadHis,UhisLoading);
Plot_Energy(UhisLoading,strainEnergyLoading,UhisAssemble,strainEnergyAssemble);
toc
Plot_DeformedHis(viewControl,newNode,newPanel,UhisLoading,UhisAssemble);




