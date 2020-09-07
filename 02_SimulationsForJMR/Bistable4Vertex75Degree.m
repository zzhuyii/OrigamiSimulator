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
z=9.1959*10^(-3);


node0=[0 0 z;
      -b 0 0;
      0 -b 0;
      b 0 0;
      0 b 0;
      -a -a -z;
      a -a -z;
      a a -z;
      -a a -z;];
  
panel0{1}=[2 1 5 9];
panel0{2}=[1 4 8 5];
panel0{3}=[6 3 1 2];
panel0{4}=[3 7 4 1];


u=node0(2,:)-node0(1,:);
v=node0(3,:)-node0(1,:);
cosTheta = dot(u,v)/(norm(u)*norm(v));
thetaInDegrees = acosd(cosTheta);
missingBeta=4*(90-thetaInDegrees)


%% Setting up the plots for display

viewControl=zeros(10,1);

% View1: View angle 1
viewControl(1)=15; 

% View2: View angle 2
viewControl(2)=30; 

% Vsize: displayed axis range 
viewControl(3)=60*10^(-3); 

% Vratio: ratio of displayed negative axis range versus the positive axis range
viewControl(4)=1; 

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
creaseWidthMat(2)=6*10^(-3);
creaseWidthMat(3)=6*10^(-3);
creaseWidthMat(5)=6*10^(-3);
creaseWidthMat(10)=6*10^(-3);
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
modelMechanicalConstant{1}=2852*10^6; 

% CreaseE: Young's modulus of creases
modelMechanicalConstant{2}=25*10^6; 

% PanelPoisson: Poisson ratio of panel
modelMechanicalConstant{3}=0.4; 

% CreasePoisson: Poisson ratio of crease
modelMechanicalConstant{4}=0.4; 

% PanelThick: thickness of panel;
% This is a vector storeing the thicknes of panels. 
modelMechanicalConstant{5}=[1;1;1;1]*2.2*10^(-3); 

% thickness of creases;
creaesThicknessMat=zeros(oldCreaseNum,1);
creaesThicknessMat(2)=1*10^(-3);
creaesThicknessMat(3)=1*10^(-3);
creaesThicknessMat(5)=1*10^(-3);
creaesThicknessMat(10)=1*10^(-3);
modelMechanicalConstant{6}=creaesThicknessMat; 


% DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.
modelMechanicalConstant{7}=4; 

% panelW
% this is an averaged creaseW. used to calculate panel bending stiffness
panelW=6*10^-3;
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



%% Input support information

% define support information
supp=[2,1,1,1;
      5,1,1,1;
      11,1,1,1;
      16,1,1,1;];
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

%% loading information for loading process  
  
loadForce=10;
load=[29,0,0,loadForce;
      30,0,0,loadForce;
      31,0,0,loadForce;
      32,0,0,loadForce;];


%% Nonlinear Solver for loading
loadConstant=zeros(4,1);

loadConstant(1)=30; % IncreStep
loadConstant(2)=10^-6; % Tor
loadConstant(3)=50; % iterMax
loadConstant(4)=1; % LambdaBar

AssembleNode=newNode;
A=size(newNode);
UhisAssemble=zeros(1,A(1),A(2));
U=zeros(size(newNode));
strainEnergyAssemble=zeros(1,4);
[sprTargetZeroStrain]=Spr_Theta(U,sprIJKL,newNode);

% if displacement controled method is used, use the following to set up the
% controlling displacement entry
dispControler=[29,3];

[U,UhisLoading,loadhis,strainEnergyLoading,...
    nodeForce,loadForce,contactForce]...
    =Solver_LoadingDC(panel0,newNode2OldNode,oldCreaseNum,...
    newNode,barConnect,barType,barArea,barLength,...
    sprIJKL,sprK,sprTargetZeroStrain,creaseRef,load,supportInfo,U,...
    loadConstant,modelGeometryConstant,modelMechanicalConstant,dispControler);
 

%% Plotting the results
deformNode=U+newNode;
Plot_DeformedShape(viewControl,AssembleNode,deformNode,newPanel);
Plot_LoadHis(loadhis,UhisLoading);
Plot_Energy(UhisLoading,strainEnergyLoading,UhisAssemble,strainEnergyAssemble);
Plot_DeformedHis(viewControl,newNode,newPanel,UhisLoading,UhisAssemble);
