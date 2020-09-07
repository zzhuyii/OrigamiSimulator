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
node0=[0 0 0;
      0 3*a 0;
      0 a*6 0;
      a*3 0 0;
      a*3 a*3 0;
      a*3 a*3.5 0;
      a*3 a*4.5 0;
      a*3 a*6 0;
      a*4.5 a*4 0;
      a*5 a*3 0;
      a*5.5 a*5 0;
      a*6 0 0;
      a*6 a*3 0;];
  
panel0{1}=[1 4 5 2];
panel0{2}=[2 5 6 7 8 3];
panel0{3}=[4 12 13 10 5];
panel0{4}=[6 9 7];
panel0{5}=[13 11 10];


%% Setting up the plots for display

viewControl=zeros(4,1);

 % View1: View angle 1
viewControl(1)=45;

 % View2: View angle 2
viewControl(2)=45;

% Vsize: displayed axis range 
viewControl(3)=150*10^(-3); 

% Vratio: ratio of displayed negative axis range versus the positive
viewControl(4)=0.2; 



%% Generate the geometry of compliant crase

% identify the crease based on the original input;
[oldCreaseNum,oldCreaseConnect,oldCreaseType]=Mesh_IdentifyCrease(node0,panel0);

% Plot the original meshing for inspection
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
creaseWidthMat(3)=8*10^(-3);
creaseWidthMat(4)=8*10^(-3);
creaseWidthMat(7)=8*10^(-3);
creaseWidthMat(12)=8*10^(-3);
modelGeometryConstant{3}=creaseWidthMat; 

% generate the geometry of system
[newNode,newPanel,barType,barConnect,...
    sprIJKL,type1BarNum,panelInnerBarStart,centerNodeStart,...
    newNode2OldNode,newCrease2OldCrease,newPanel2OldPanel,newPanelNum] ...
    =Mesh_CompliantCreaseGeometry(node0,panel0,...
    oldCreaseNum,oldCreaseConnect,oldCreaseType,...
    modelGeometryConstant);

% We need to manually adjust some nodal coordinates for this case
creaseW=8*10^(-3);
newNode(15,2)=newNode(15,2)-0.5*creaseW;
newNode(14,2)=newNode(14,2)-0.5*creaseW;    
newNode(32,2)=newNode(32,2)-0*creaseW;
newNode(33,2)=newNode(33,2)+0.125*creaseW;
newNode(31,2)=newNode(31,2)+0.25*creaseW;
newNode(21,2)=newNode(21,2)+0.5*creaseW;
newNode(19,2)=newNode(19,2)+0.5*creaseW;    
newNode(23,2)=newNode(23,2)-0.25*creaseW;
newNode(24,2)=newNode(24,2)-0.125*creaseW;    


newNode(6,1)=newNode(6,1)-0.5*creaseW;
newNode(25,1)=newNode(25,1)-0.25*creaseW;
newNode(27,1)=newNode(27,1)-0.125*creaseW;    
newNode(7,1)=newNode(7,1)-0.5*creaseW;
newNode(8,1)=newNode(8,1)-0.5*creaseW;
newNode(9,1)=newNode(9,1)-0.5*creaseW; 

newNode(28,1)=newNode(28,1)+0.125*creaseW;
newNode(29,1)=newNode(29,1)+0.125*creaseW;
newNode(30,1)=newNode(30,1)+0.125*creaseW;  
newNode(16,1)=newNode(16,1)+0.5*creaseW;
newNode(18,1)=newNode(18,1)+0.5*creaseW;


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




%% Assign mechanical properties

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
modelMechanicalConstant{5}=[1;1;1;1;1]*200*10^(-6); 

% thickness of creases;
creaesThicknessMat=zeros(oldCreaseNum,1);
creaesThicknessMat(3)=100*10^(-6);
creaesThicknessMat(4)=100*10^(-6);
creaesThicknessMat(7)=100*10^(-6);
creaesThicknessMat(12)=100*10^(-6);
modelMechanicalConstant{6}=creaesThicknessMat; 

% DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.
modelMechanicalConstant{7}=4; 

% panelW
% this is an averaged creaseW. used to calculate panel bending stiffness
panelW=8*10^-3;
modelMechanicalConstant{8}=panelW;




%% panel contact related input

% contact open
% 1: means consider panel contact;
% 0: means ignore panel contact
modelMechanicalConstant{9}= 1; 

% ke: used to scale the magnitude of potentil
modelMechanicalConstant{10}=0.02; 

% d0edge: d0 for points at the edge
modelMechanicalConstant{11}=4*10^(-3);

% d0center: d0 for points at the center
modelMechanicalConstant{12}=2*10^(-3);




%% Assign zero strain position for creases during self-folding

% 0-2pi, This matrix can be used to manually set the zero energy rotation
% angle of the crease hinge
rotationZeroStrain=pi*ones(oldCreaseNum,1);

% Folding Sequence indicate which crease will be folded first
foldingSequence=ones(oldCreaseNum,1);

% Target zero strain folding angle of creases;
ratio=0.5;
rotationZeroStrain(3)=pi+ratio*pi;
rotationZeroStrain(4)=pi+ratio*pi;
foldingSequence(3)=1;
foldingSequence(4)=1;

ratio2=0.5;
rotationZeroStrain(7)=pi+ratio2*pi;
rotationZeroStrain(12)=pi+ratio2*pi;
foldingSequence(7)=2;
foldingSequence(12)=2;

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
      3,1,1,1;
      4,1,1,1;];
supportInfo{1}=supp;    

% NonRigidSupport
% 0 means the non rigid support is not activated.
% 1 means the non rigid support is activated.
supportInfo{2}=0; 

% first column stores node number
% second column stores direction
% third column stores stiffness
ksup=500;
suppElastic=[1,3,ksup;];
supportInfo{3}=suppElastic;



%% information for Self folding process
loadForce=0;
load=[12,0,0,-loadForce;];
  
  

%% Nonlinear Solver for Assemble

assembleConstant=zeros(3,1);

% increStep
assembleConstant(1)=50; 

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

loadForce=0.5*10^(-3);
load=[12,loadForce,0,0;];


%% Nonlinear Solver for loading

loadConstant=zeros(4,1);

% increStep
loadConstant(1)=50; 

% tor
loadConstant(2)=10^-8; 

% iterMax
loadConstant(3)=50; 

% lambdaBar
loadConstant(4)=1; 

[U,UhisLoading,loadhis,strainEnergyLoading,nodeForce,loadForce,lockForce]...
    =Solver_LoadingNR(panel0,newNode2OldNode,oldCreaseNum,...
    newNode,barConnect,barType,barArea,barLength, ...
    sprIJKL,sprK,sprTargetZeroStrain,creaseRef,load,supportInfo, ...
    U,loadConstant,modelGeometryConstant,modelMechanicalConstant);



%% Plotting the results
deformNode=U+newNode;
Plot_DeformedHis(viewControl,newNode,newPanel,UhisLoading,UhisAssemble);


