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
node0=[0 0 0;
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
      
panel0{1}=[1 2 3];
panel0{2}=[1 3 4];
panel0{3}=[1 4 5];
panel0{4}=[1 5 6];
panel0{5}=[1 6 7];
panel0{6}=[1 7 8];
panel0{7}=[1 8 9];
panel0{8}=[1 9 2];

panel0{9}=[2 10 11 3];
panel0{10}=[3 12 13 4];
panel0{11}=[4 14 15 5];
panel0{12}=[5 16 17 6];
panel0{13}=[6 18 19 7];
panel0{14}=[7 20 21 8];
panel0{15}=[8 22 23 9];
panel0{16}=[9 24 25 2];

panel0{17}=[2 26 10];
panel0{18}=[3 11 27];
panel0{19}=[3 27 12];
panel0{20}=[4 13 28];
panel0{21}=[4 28 14];
panel0{22}=[5 15 29];
panel0{23}=[5 29 16];
panel0{24}=[6 17 30];
panel0{25}=[6 30 18];
panel0{26}=[7 19 31];
panel0{27}=[7 31 20];
panel0{28}=[8 21 32];
panel0{29}=[8 32 22];
panel0{30}=[9 23 33];
panel0{31}=[9 33 24];
panel0{32}=[2 25 26];


panel0{33}=[10 34 11];
panel0{34}=[12 35 13];
panel0{35}=[14 36 15];
panel0{36}=[16 37 17];
panel0{37}=[18 38 19];
panel0{38}=[20 39 21];
panel0{39}=[22 40 23];
panel0{40}=[24 41 25];


%% Setting up the plots for display

viewControl=zeros(10,1);

% View1: View angle 1
viewControl(1)=300; 

% View2: View angle 2
viewControl(2)=30; 

% Vsize: displayed axis range 
viewControl(3)=50*10^(-3); 

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

creaseWidthMat(1)=1.5*10^(-3);
creaseWidthMat(2)=1.5*10^(-3);
creaseWidthMat(14)=1.5*10^(-3);
creaseWidthMat(12)=1.5*10^(-3);
creaseWidthMat(10)=1.5*10^(-3);
creaseWidthMat(8)=1.5*10^(-3);
creaseWidthMat(6)=1.5*10^(-3);
creaseWidthMat(4)=1.5*10^(-3);

creaseWidthMat(3)=1.5*10^(-3);
creaseWidthMat(5)=1.5*10^(-3);
creaseWidthMat(7)=1.5*10^(-3);
creaseWidthMat(9)=1.5*10^(-3);
creaseWidthMat(11)=1.5*10^(-3);
creaseWidthMat(13)=1.5*10^(-3);
creaseWidthMat(15)=1.5*10^(-3);
creaseWidthMat(16)=1.5*10^(-3);

creaseWidthMat(41)=1.5*10^(-3);
creaseWidthMat(43)=1.5*10^(-3);
creaseWidthMat(46)=1.5*10^(-3);
creaseWidthMat(49)=1.5*10^(-3);
creaseWidthMat(52)=1.5*10^(-3);
creaseWidthMat(55)=1.5*10^(-3);
creaseWidthMat(58)=1.5*10^(-3);
creaseWidthMat(61)=1.5*10^(-3);

creaseWidthMat(40)=1.5*10^(-3);
creaseWidthMat(17)=1.5*10^(-3);
creaseWidthMat(19)=1.5*10^(-3);
creaseWidthMat(20)=1.5*10^(-3);
creaseWidthMat(22)=1.5*10^(-3);
creaseWidthMat(23)=1.5*10^(-3);
creaseWidthMat(25)=1.5*10^(-3);
creaseWidthMat(26)=1.5*10^(-3);
creaseWidthMat(28)=1.5*10^(-3);
creaseWidthMat(29)=1.5*10^(-3);
creaseWidthMat(31)=1.5*10^(-3);
creaseWidthMat(32)=1.5*10^(-3);
creaseWidthMat(34)=1.5*10^(-3);
creaseWidthMat(35)=1.5*10^(-3);
creaseWidthMat(37)=1.5*10^(-3);
creaseWidthMat(38)=1.5*10^(-3);

creaseWidthMat(18)=1.5*10^(-3);
creaseWidthMat(21)=1.5*10^(-3);
creaseWidthMat(24)=1.5*10^(-3);
creaseWidthMat(27)=1.5*10^(-3);
creaseWidthMat(30)=1.5*10^(-3);
creaseWidthMat(33)=1.5*10^(-3);
creaseWidthMat(36)=1.5*10^(-3);
creaseWidthMat(39)=1.5*10^(-3);

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
modelMechanicalConstant{2}=20*10^6; 

% PanelPoisson: Poisson ratio of panel
modelMechanicalConstant{3}=0.3; 

% CreasePoisson: Poisson ratio of crease
modelMechanicalConstant{4}=0.3; 

% PanelThick: thickness of panel;
% This is a vector storeing the thicknes of panels. 
modelMechanicalConstant{5}=ones(40,1)*1000*10^(-6); 

% thickness of creases;
% crease width to generate the topology of the compliant 
% crease origami;
creaesThicknessMat=zeros(oldCreaseNum,1);

creaesThicknessMat(1)=300*10^(-6);
creaesThicknessMat(2)=300*10^(-6);
creaesThicknessMat(14)=300*10^(-6);
creaesThicknessMat(12)=300*10^(-6);
creaesThicknessMat(10)=300*10^(-6);
creaesThicknessMat(8)=300*10^(-6);
creaesThicknessMat(6)=300*10^(-6);
creaesThicknessMat(4)=300*10^(-6);

creaesThicknessMat(3)=300*10^(-6);
creaesThicknessMat(5)=300*10^(-6);
creaesThicknessMat(7)=300*10^(-6);
creaesThicknessMat(9)=300*10^(-6);
creaesThicknessMat(11)=300*10^(-6);
creaesThicknessMat(13)=300*10^(-6);
creaesThicknessMat(15)=300*10^(-6);
creaesThicknessMat(16)=300*10^(-6);

creaesThicknessMat(41)=300*10^(-6);
creaesThicknessMat(43)=300*10^(-6);
creaesThicknessMat(46)=300*10^(-6);
creaesThicknessMat(49)=300*10^(-6);
creaesThicknessMat(52)=300*10^(-6);
creaesThicknessMat(55)=300*10^(-6);
creaesThicknessMat(58)=300*10^(-6);
creaesThicknessMat(61)=300*10^(-6);

creaesThicknessMat(40)=300*10^(-6);
creaesThicknessMat(17)=300*10^(-6);
creaesThicknessMat(19)=300*10^(-6);
creaesThicknessMat(20)=300*10^(-6);
creaesThicknessMat(22)=300*10^(-6);
creaesThicknessMat(23)=300*10^(-6);
creaesThicknessMat(25)=300*10^(-6);
creaesThicknessMat(26)=300*10^(-6);
creaesThicknessMat(28)=300*10^(-6);
creaesThicknessMat(29)=300*10^(-6);
creaesThicknessMat(31)=300*10^(-6);
creaesThicknessMat(32)=300*10^(-6);
creaesThicknessMat(34)=300*10^(-6);
creaesThicknessMat(35)=300*10^(-6);
creaesThicknessMat(37)=300*10^(-6);
creaesThicknessMat(38)=300*10^(-6);

creaesThicknessMat(18)=300*10^(-6);
creaesThicknessMat(21)=300*10^(-6);
creaesThicknessMat(24)=300*10^(-6);
creaesThicknessMat(27)=300*10^(-6);
creaesThicknessMat(30)=300*10^(-6);
creaesThicknessMat(33)=300*10^(-6);
creaesThicknessMat(36)=300*10^(-6);
creaesThicknessMat(39)=300*10^(-6);
modelMechanicalConstant{6}=creaesThicknessMat; 


% DiagonalRate:
% Diagonal Rate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.
modelMechanicalConstant{7}=4; 

% panelW
% this is an averaged creaseW. used to calculate panel bending stiffness
panelW=1.5*10^-3;
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

GlobalRate=1;

ratio=-0.4*GlobalRate;
rotationZeroStrain(3)=pi-ratio*pi;
rotationZeroStrain(5)=pi-ratio*pi;
rotationZeroStrain(7)=pi-ratio*pi;
rotationZeroStrain(9)=pi-ratio*pi;
rotationZeroStrain(11)=pi-ratio*pi;
rotationZeroStrain(13)=pi-ratio*pi;
rotationZeroStrain(15)=pi-ratio*pi;
rotationZeroStrain(16)=pi-ratio*pi;

ratio1=-0.4*GlobalRate;
rotationZeroStrain(41)=pi+ratio1*pi;
rotationZeroStrain(43)=pi-ratio1*pi;
rotationZeroStrain(46)=pi-ratio1*pi;
rotationZeroStrain(49)=pi-ratio1*pi;
rotationZeroStrain(52)=pi-ratio1*pi;
rotationZeroStrain(55)=pi-ratio1*pi;
rotationZeroStrain(58)=pi-ratio1*pi;
rotationZeroStrain(61)=pi-ratio1*pi;

rotationZeroStrain(40)=pi-ratio1*pi;
rotationZeroStrain(17)=pi-ratio1*pi;
rotationZeroStrain(19)=pi-ratio1*pi;
rotationZeroStrain(20)=pi-ratio1*pi;
rotationZeroStrain(22)=pi-ratio1*pi;
rotationZeroStrain(23)=pi-ratio1*pi;
rotationZeroStrain(25)=pi-ratio1*pi;
rotationZeroStrain(26)=pi-ratio1*pi;
rotationZeroStrain(28)=pi-ratio1*pi;
rotationZeroStrain(29)=pi-ratio1*pi;
rotationZeroStrain(31)=pi-ratio1*pi;
rotationZeroStrain(32)=pi-ratio1*pi;
rotationZeroStrain(34)=pi-ratio1*pi;
rotationZeroStrain(35)=pi-ratio1*pi;
rotationZeroStrain(37)=pi-ratio1*pi;
rotationZeroStrain(38)=pi-ratio1*pi;

ratio2=-0.4*GlobalRate;
rotationZeroStrain(18)=pi-ratio2*pi;
rotationZeroStrain(21)=pi-ratio2*pi;
rotationZeroStrain(24)=pi-ratio2*pi;
rotationZeroStrain(27)=pi-ratio2*pi;
rotationZeroStrain(30)=pi-ratio2*pi;
rotationZeroStrain(33)=pi-ratio2*pi;
rotationZeroStrain(36)=pi-ratio2*pi;
rotationZeroStrain(39)=pi-ratio2*pi;

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
      3,1,1,1;];
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

loadForce=0*10^(-3);
load=[1,0,0,-loadForce;];

%% Nonlinear Solver for Assemble
assembleConstant=zeros(3,1);
assembleConstant(1)=80; % IncreStep
assembleConstant(2)=5*10^-5; % Tor
assembleConstant(3)=50; % iterMax

[U,UhisAssemble,strainEnergyAssemble]=Solver_Assemble(...
    panel0,newNode2OldNode,newNode,barConnect,barType,barLength,...
    barArea,sprIJKL,sprK,sprTargetZeroStrain,sprFoldingSequence,...
    creaseRef,oldCreaseNum,assembleConstant,...
    modelGeometryConstant,modelMechanicalConstant,...
    supportInfo,load);

assembleNode=U+newNode;
Plot_DeformedShapeOnly(viewControl,assembleNode,newPanel)


%% Support and loading information for loading process
loadForce=0.1*10^(-3);
load=[1,0,0,-loadForce;];


%% Nonlinear Solver for loading
loadConstant=zeros(4,1);
loadConstant(1)=1; % IncreStep
loadConstant(2)=10^-8; % Tor
loadConstant(3)=50; % iterMax
loadConstant(4)=0.01; % LambdaBar


[U,UhisLoading,loadhis,strainEnergyLoading,nodeForce,loadForce,lockForce]...
    =Solver_LoadingNR(panel0,newNode2OldNode,oldCreaseNum,...
    newNode,barConnect,barType,barArea,barLength, ...
    sprIJKL,sprK,sprTargetZeroStrain,creaseRef,load,supportInfo, ...
    U,loadConstant,modelGeometryConstant,modelMechanicalConstant);



%% Plotting the results
deformNode=U+newNode;
Plot_Energy(UhisLoading,strainEnergyLoading,UhisAssemble,strainEnergyAssemble);

toc

angleNum=zeros(4,1);
angleNum(1)=543;
angleNum(2)=544;
% Midline Crease Num
angleNum(3)=60;
angleNum(4)=64;
% Crease Num at the two panels
[angleHis]=Plot_DeformedHisAndCalcTheta(viewControl,newNode,newPanel,...
    UhisLoading,UhisAssemble,sprIJKL,angleNum);

