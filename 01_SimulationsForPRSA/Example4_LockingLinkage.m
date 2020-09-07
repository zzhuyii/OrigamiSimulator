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
node0=[b 0 0;
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
  
panel0{1}=[1 2 8 7];
panel0{2}=[2 3 9 8];
panel0{3}=[3 4 10 9];
panel0{4}=[4 5 11 10];
panel0{5}=[5 6 12 11];
panel0{6}=[7 8 14 13];
panel0{7}=[8 9 15 14];
panel0{8}=[9 10 16 15];
panel0{9}=[10 11 17 16];
panel0{10}=[11 12 18 17];

panel0{11}=[1 2 8 7]+[18 18 18 18];
panel0{12}=[2 3 9 8]+[18 18 18 18];
panel0{13}=[3 4 10 9]+[18 18 18 18];
panel0{14}=[4 5 11 10]+[18 18 18 18];
panel0{15}=[5 6 12 11]+[18 18 18 18];
panel0{16}=[7 8 14 13]+[18 18 18 18];
panel0{17}=[8 9 15 14]+[18 18 18 18];
panel0{18}=[9 10 16 15]+[18 18 18 18];
panel0{19}=[10 11 17 16]+[18 18 18 18];
panel0{20}=[11 12 18 17]+[18 18 18 18];


%% Setting up the plots for display

viewControl=zeros(10,1);

% View1: View angle 1
viewControl(1)=15;

% View2: View angle 2
viewControl(2)=15; 

% Vsize: displayed axis range 
viewControl(3)=70*10^(-3);

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
modelGeometryConstant{1}=3; 

% 1 means include compliant crease model
% 0 means using concentrated hinge model
modelGeometryConstant{2}=0; 

% crease width to generate the topology of the compliant 
% crease origami;
creaseWidthMat=zeros(oldCreaseNum,1);

creaseWidthMat(4)=3*10^(-3);
creaseWidthMat(7)=3*10^(-3);
creaseWidthMat(10)=3*10^(-3);
creaseWidthMat(13)=3*10^(-3);
creaseWidthMat(16)=3*10^(-3);

creaseWidthMat(3)=3*10^(-3);
creaseWidthMat(18)=3*10^(-3);
creaseWidthMat(6)=3*10^(-3);
creaseWidthMat(20)=3*10^(-3);
creaseWidthMat(9)=3*10^(-3);
creaseWidthMat(22)=3*10^(-3);
creaseWidthMat(12)=3*10^(-3);
creaseWidthMat(24)=3*10^(-3);

creaseWidthMat(4+27)=3*10^(-3);
creaseWidthMat(7+27)=3*10^(-3);
creaseWidthMat(10+27)=3*10^(-3);
creaseWidthMat(13+27)=3*10^(-3);
creaseWidthMat(16+27)=3*10^(-3);

creaseWidthMat(3+27)=3*10^(-3);
creaseWidthMat(18+27)=3*10^(-3);
creaseWidthMat(6+27)=3*10^(-3);
creaseWidthMat(20+27)=3*10^(-3);
creaseWidthMat(9+27)=3*10^(-3);
creaseWidthMat(22+27)=3*10^(-3);
creaseWidthMat(12+27)=3*10^(-3);
creaseWidthMat(24+27)=3*10^(-3);

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
modelMechanicalConstant{1}=2*10^9; 

% CreaseE: Young's modulus of creases
modelMechanicalConstant{2}=2*10^9; 

% PanelPoisson: Poisson ratio of panel
modelMechanicalConstant{3}=0.3; 

% CreasePoisson: Poisson ratio of crease
modelMechanicalConstant{4}=0.3; 

% PanelThick: thickness of panel;
% This is a vector storeing the thicknes of panels. 
modelMechanicalConstant{5}=ones(20,1)*500*10^(-6);

% thickness of creases;
creaesThicknessMat=zeros(oldCreaseNum,1);

creaesThicknessMat(4)=100*10^(-6);
creaesThicknessMat(7)=100*10^(-6);
creaesThicknessMat(10)=100*10^(-6);
creaesThicknessMat(13)=100*10^(-6);
creaesThicknessMat(16)=100*10^(-6);

creaesThicknessMat(3)=100*10^(-6);
creaesThicknessMat(18)=100*10^(-6);
creaesThicknessMat(6)=100*10^(-6);
creaesThicknessMat(20)=100*10^(-6);
creaesThicknessMat(9)=100*10^(-6);
creaesThicknessMat(22)=100*10^(-6);
creaesThicknessMat(12)=100*10^(-6);
creaesThicknessMat(24)=100*10^(-6);

creaesThicknessMat(4+27)=100*10^(-6);
creaesThicknessMat(7+27)=100*10^(-6);
creaesThicknessMat(10+27)=100*10^(-6);
creaesThicknessMat(13+27)=100*10^(-6);
creaesThicknessMat(16+27)=100*10^(-6);

creaesThicknessMat(3+27)=100*10^(-6);
creaesThicknessMat(18+27)=100*10^(-6);
creaesThicknessMat(6+27)=100*10^(-6);
creaesThicknessMat(20+27)=100*10^(-6);
creaesThicknessMat(9+27)=100*10^(-6);
creaesThicknessMat(22+27)=100*10^(-6);
creaesThicknessMat(12+27)=100*10^(-6);
creaesThicknessMat(24+27)=100*10^(-6);

modelMechanicalConstant{6}=creaesThicknessMat; 


% DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.
modelMechanicalConstant{7}=4; 

% panelW
% this is an averaged creaseW. used to calculate panel bending stiffness
panelW=3*10^-3;
modelMechanicalConstant{8}=panelW;



%% panel contact related input

% contact open
% 1: means consider panel contact;
% 0: means ignore panel contact
modelMechanicalConstant{9}= 0; 

% ke: used to scale the magnitude of potentil
modelMechanicalConstant{10}=0.002; 

% d0edge: d0 for points at the edge
modelMechanicalConstant{11}=1*10^-3; 

% d0center: d0 for points at the center
modelMechanicalConstant{12}=1*10^-3;  



%% Assign zero strain position for creases during self-folding

% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
rotationZeroStrain=pi*ones(oldCreaseNum,1);

foldingSequence=ones(oldCreaseNum,1);
% Folding Sequence indicate which crease will be folded first

ratio=0.5;
rotationZeroStrain(4)=pi+ratio*pi;
rotationZeroStrain(7)=pi-ratio*pi;
rotationZeroStrain(10)=pi+ratio*pi;
rotationZeroStrain(13)=pi-ratio*pi;
rotationZeroStrain(16)=pi+ratio*pi;

rotationZeroStrain(3)=pi-ratio*pi;
rotationZeroStrain(18)=pi-ratio*pi;
rotationZeroStrain(6)=pi+ratio*pi;
rotationZeroStrain(20)=pi+ratio*pi;
rotationZeroStrain(9)=pi+ratio*pi;
rotationZeroStrain(22)=pi+ratio*pi;
rotationZeroStrain(12)=pi-ratio*pi;
rotationZeroStrain(24)=pi-ratio*pi;

rotationZeroStrain(4+27)=pi+ratio*pi;
rotationZeroStrain(7+27)=pi-ratio*pi;
rotationZeroStrain(10+27)=pi+ratio*pi;
rotationZeroStrain(13+27)=pi-ratio*pi;
rotationZeroStrain(16+27)=pi+ratio*pi;

rotationZeroStrain(3+27)=pi-ratio*pi;
rotationZeroStrain(18+27)=pi-ratio*pi;
rotationZeroStrain(6+27)=pi+ratio*pi;
rotationZeroStrain(20+27)=pi+ratio*pi;
rotationZeroStrain(9+27)=pi+ratio*pi;
rotationZeroStrain(22+27)=pi+ratio*pi;
rotationZeroStrain(12+27)=pi-ratio*pi;
rotationZeroStrain(24+27)=pi-ratio*pi;

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



%% Input information of support for Assemble

% define support information
supp=[1,1,1,1;
      2,1,1,1;
      13,0,0,1;
      14,0,0,1;
      19,1,1,1;
      20,1,1,1;
      31,0,0,1;
      32,0,0,1;];
supportInfo{1}=supp;    

  
% NonRigidSupport
% 0 means the non rigid support is not activated.
% 1 means the non rigid support is activated.
supportInfo{2}=0; 


% first column stores node number
% second column stores direction
% third column stores stiffness
ksup=500;
suppElastic=[19,3,ksup;
             20,3,ksup;
             31,3,ksup;
             32,3,ksup;
             23,3,ksup;
             24,3,ksup;
             35,3,ksup;
             36,3,ksup;];
supportInfo{3}=suppElastic;



%% information for Self folding process
loadForce=0;
load=[10,0,0,-loadForce;];
  
%% Nonlinear Solver for Assemble

assembleConstant=zeros(3,1);

% increStep
assembleConstant(1)=30;

% tor
assembleConstant(2)=10^-6; 

% iterMax
assembleConstant(3)=50; 

[U,UhisAssemble,strainEnergyAssemble]=Solver_Assemble(...
    panel0,newNode2OldNode,newNode,barConnect,barType,barLength,...
    barArea,sprIJKL,sprK,sprTargetZeroStrain,sprFoldingSequence,...
    creaseRef,oldCreaseNum,assembleConstant,...
    modelGeometryConstant,modelMechanicalConstant,...
    supportInfo,load);

assembleNode=U+newNode;
Plot_DeformedShapeOnly(viewControl,assembleNode,newPanel)

%% Support and loading information for loading process

supp=[2,1,1,1;
      14,1,1,1;
      5,1,1,1;
      17,1,1,1;
      32,1,1,0;
      20,1,1,0;
      23,1,1,0;
      35,1,1,0;];
supportInfo{1}=supp;   

% NonRigidSupport
% 0 means the non rigid support is not activated.
% 1 means the non rigid support is activated.
supportInfo{2}=1; 


% first column stores node number
% second column stores direction
% third column stores stiffness
ksup=5;
suppElastic=[19,3,ksup;
             20,3,ksup;
             31,3,ksup;
             32,3,ksup;
             23,3,ksup;
             24,3,ksup;
             35,3,ksup;
             36,3,ksup;];
supportInfo{3}=suppElastic;

%% loading information

loadForce=5*10^(-3);
load=[21,0,0,loadForce;
      22,0,0,loadForce;
      33,0,0,loadForce;
      34,0,0,loadForce;]; 
  
% contact open
% 1: means consider panel contact;
% 0: means ignore panel contact
modelMechanicalConstant{9}= 1; 

%% Nonlinear Solver for loading

loadConstant=zeros(4,1);

% increStep
loadConstant(1)=30; 

% tor
loadConstant(2)=10^-6; 

% iterMax
loadConstant(3)=50; 

% lambdaBar
loadConstant(4)=1; 


[U,UhisLoading,loadHis,strainEnergyLoading,nodeForce,loadForce,lockForce]...
    =Solver_LoadingNR(panel0,newNode2OldNode,oldCreaseNum,...
    newNode,barConnect,barType,barArea,barLength, ...
    sprIJKL,sprK,sprTargetZeroStrain,creaseRef,load,supportInfo, ...
    U,loadConstant,modelGeometryConstant,modelMechanicalConstant);


%% Plotting the results
deformNode=U+newNode;
Plot_DeformedHis(viewControl,newNode,newPanel,UhisLoading,UhisAssemble);
