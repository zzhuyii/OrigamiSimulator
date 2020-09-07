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

a=1*10^(-3);
b=100*10^-6;
node0=[0 0 0;
      a 0 0;
      2*a 0 0;
      0 a 0;
      a a 0;
      2*a a 0;
      0*a,-0.6*a,-b;
      0*a,1.6*a,-b;
      2.5*a,-0.6*a,-b;
      2.5*a,1.6*a,-b;];
  
panel0{1}=[1 2 5 4];
panel0{2}=[2 3 6 5];
panel0{3}=[7 8 10 9];

%% Setting up the plots for display

viewControl=zeros(10,1);

% View1: View angle 1
viewControl(1)=45; 

% View2: View angle 2
viewControl(2)=45; 

% Vsize: displayed axis range 
viewControl(3)=3*10^(-3); 

% Vratio: ratio of displayed negative axis range versus the positive axis range
viewControl(4)=0.2; 



%% Generate the geometry of compliant crease

% Here we identify the creases based on the original input;
[oldCreaseNum,oldCreaseConnect,oldCreaseType]=Mesh_IdentifyCrease(node0,panel0);

% Plot the original meshing for inspection;
Plot_OriginalMeshing(node0,panel0,oldCreaseNum,oldCreaseConnect,viewControl)

% Flag2D3D is used to determine how the additional crease structuer is
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
creaseWidthMat(3)=400*10^(-6); 
modelGeometryConstant{3}=creaseWidthMat; 

% generate the geometry of system
[newNode,newPanel,barType,barConnect,...
    sprIJKL,type1BarNum,panelInerBarStart,centerNodeStart,...
    newNode2OldNode,newCrease2OldCrease,newPanel2OldPanel,panelNum] ...
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
modelGeometryConstant{4}=panelInerBarStart; 
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
modelMechanicalConstant{5}=[500;20;500]*10^(-6); 

% thickness of creases;
creaesThicknessMat=zeros(oldCreaseNum,1);
creaesThicknessMat(3)=1*10^(-6);
modelMechanicalConstant{6}=creaesThicknessMat;

% diagonalRate:
% DiagonalRate is the factor that determine how much diagonal springs are
% stiffer than horizontal ones.
modelMechanicalConstant{7}=50; 

% panelW
% this is an averaged creaseW. used to calculate panel bending stiffness
panelW=creaseWidthMat(3);
modelMechanicalConstant{8}=panelW;




%% panel contact related input

% contact open
% 1: means consider panel contact;
% 0: means ignore panel contact
modelMechanicalConstant{9}=1;

% ke: used to scale the magnitude of potentil
modelMechanicalConstant{10}=0.0001; 

 % d0edge: d0 for points at the edge
modelMechanicalConstant{11}=20*(10^(-6));

% d0center: d0 for points at the center
modelMechanicalConstant{12}=20*(10^(-6));



%% Assign zero strain position for creases during self-folding

% 0-2pi, This matrix can be used to manually set the target zero energy 
% angle of the crease hinge
rotationZeroStrain=pi*ones(oldCreaseNum,1);

% Folding Sequence indicate which crease will be folded first
foldingSequence=ones(oldCreaseNum,1);

% Target zero strain folding angle of creases;
ratio=0;
rotationZeroStrain(3)=pi+ratio*pi;

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
supp=[1,0,0,1;
      2,0,0,1;
      3,0,1,1;
      4,1,1,1;
      9,1,1,1;
      10,1,1,1;
      11,1,1,1;
      12,1,1,1;];
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
% This is used to include the effects of gravity

loadConstant=zeros(4,1);

% increStep
loadConstant(1)=80; 

% tor
loadConstant(2)=10^-6; 

% iterMax
loadConstant(3)=50; 

% lambdaBar
loadConstant(4)=0.01; 



% Here we divide the magnitude by the step number so that after the newton
% method we get the desired gravity we want.
loadMag=(40*10^(-9))/loadConstant(1);
load=[5,0,0,-loadMag;
      6,0,0,-loadMag;
      7,0,0,-loadMag;
      8,0,0,-loadMag;
      17,0,0,-2*loadMag;];

U=zeros(size(newNode));

% solve the geometry under gravity
[U,UhisLoading,loadHis,StrainEnergyLoading,NodeForce,LoadForce,lockForce]...
    =Solver_LoadingNR(panel0,newNode2OldNode,oldCreaseNum,...
    newNode,barConnect,barType,barArea,barLength, ...
    sprIJKL,sprK,sprTargetZeroStrain,creaseRef,load,supportInfo,U,...
    loadConstant,modelGeometryConstant,modelMechanicalConstant);

gravityShape=newNode+U;
Plot_DeformedShape(viewControl,newNode,gravityShape,newPanel);

% Get the load vector with the original magnitude for latter use
loadMag=(40*10^(-9));
load=[5,0,0,-loadMag;
      6,0,0,-loadMag;
      7,0,0,-loadMag;
      8,0,0,-loadMag;
      17,0,0,-2*loadMag;];
  
  
%% Generate Thermal Conductivity Matrix and perform analysis

% Thermal Conductivity of panel 
modelThermalConstant{1} = [1.3;0.3;1.3]; 

% Thermal Conductivity of crease
modelThermalConstant{2} = 0.3;

% Thermal Conductivity of submerged environment
modelThermalConstant{3} = 0.026; 

% Thickness of the submerged environment at RT for panels
% This value is allowed to be different to consider the anchorage.
modelThermalConstant{4} = [1;1;1]*1000*10^(-6); 

% Thickness of the submerged environment at RT for creases
modelThermalConstant{5} = 1000*10^(-6); 

% Number of Air layers used
modelThermalConstant{6} = 10; 

% Thermal dissipation angle
modelThermalConstant{7} = 20/180*3.14; 

% Temperature of the submerged environment
modelThermalConstant{8} = 21; 


% Define the BC for thermal conduction, to confine the maximum air
% thickness. The vector stores the number of panels that serves as the BC
% for heat transfer (ie. those that are Si wafers).
thermalBoundaryPanels=[3];
modelThermalConstant{9}=thermalBoundaryPanels;

% This vector stores extra nodes that should be at RT to adjust the BC.
roomTempNode=[1;4];
modelThermalConstant{10}=roomTempNode;



% differential thermal expansion coefficients
modelTimoshenkoConstant{1} = (52-14)*10^(-6);  

% Young's modulus of bimorph material 1
modelTimoshenkoConstant{2}= 39.5*10^9; 

% Young's modulus of bimorph material 2
modelTimoshenkoConstant{3} = 2*10^9;

% thickness of bimorph material 1
modelTimoshenkoConstant{4}= 0.2*10^-6;

% thickness of bimorph material 2
modelTimoshenkoConstant{5}= 0.8*10^-6;




% Assemble the thermal conductivity matrix
[thermalMat,thermalNodeNum]=Thermal_AssembleConductMat...
    (newNode,newPanel,barConnect,barLength,barType,...
    U,modelThermalConstant,modelMechanicalConstant,...
    newCrease2OldCrease,newPanel2OldPanel);



% define input vector of energy
qtotal=18/1000; % total input energy with unit W
qin=zeros((modelThermalConstant{6}+1)*thermalNodeNum,1);
qin(2)=1/6*1/2*qtotal;
qin(3)=1/6*1/2*qtotal;
qin(5)=1/6*1/2*qtotal;
qin(8)=1/6*1/2*qtotal;
qin(13)=2/3*1/6*qtotal;
qin(14)=2/3*1/6*qtotal;
qin(15)=2/3*2/3*qtotal;



% Define the analyses parameter for electro-thermal induced mehcnaical
% behavior of the origami

assembleConstant=zeros(3,1);

% mechanical increStep in thermal step
assembleConstant(1)=1; 

% tor
assembleConstant(2)=5*10^-7; 

% iterMax
assembleConstant(3)=25; 


% number of steps for solving the thermal conductivity
thermalStep=80; 

% set up storage matrix
temperatureHistory=zeros(thermalNodeNum,thermalStep);
UhisThermal=zeros(thermalStep,thermalNodeNum,3);
energyHisThermal=zeros(thermalStep,4);

tempNode=gravityShape;
qhis=zeros(thermalStep,1);
foldHis=zeros(thermalStep,1);
tempHis=zeros(thermalStep,1);
for i=1:thermalStep
    
    % linearly increment of input energy
    qtemp=i/thermalStep*qin;
    qhis(i)=i/thermalStep*qtotal;
    A=size(roomTempNode);
    Nrt=A(1);
        
    % Solve for the temperature profile
    [T,indexArray]=Thermal_SolveTemperature(qtemp,thermalMat,thermalNodeNum,...
        modelThermalConstant);
    temperatureHistory(indexArray,i)=T(1:(thermalNodeNum-Nrt));
    temperatureHistory(roomTempNode,i)=modelThermalConstant{8}*ones(Nrt,1);
    T=temperatureHistory(:,i);
    
    % Update the "SprTargetZeroStrain" with Timoshenko Model
    Tave=2/3*1/3*(T(13)+T(14)+T(15))+...
        1/3*1/4*(T(2)+T(3)+T(5)+T(8));
    rot=Thermal_Timoshenko(Tave,modelTimoshenkoConstant,...
        modelThermalConstant,creaseWidthMat(3));
    tempHis(i)=Tave;
    
    % Initialize the vectors again
    % 0-2pi, This matrix can be used to manually set the target zero energy 
    % angle of the crease hinge
    rotationZeroStrain=pi*ones(oldCreaseNum,1);

    % Folding Sequence indicate which crease will be folded first
    foldingSequence=ones(oldCreaseNum,1);
    
    % Maximum number of loop needed for sequantial folding
    totalFoldingNum=max(foldingSequence);
        
    rotationZeroStrain(3)=pi+rot;
    modelMechanicalConstant{13}=rotationZeroStrain;

    [barArea,sprK,sprTargetZeroStrain,sprFoldingSequence]...
        =Mesh_MechanicalProperty(modelMechanicalConstant,...
        modelGeometryConstant,oldCreaseType,...
        oldCreaseNum,creaseRef,barLength,panel0,...
        barConnect,newNode);
   
    [U,Uhis,strainEnergyAssemble]=Solver_Assemble(panel0,newNode2OldNode,...
        tempNode,barConnect,barType,barLength,barArea,...
        sprIJKL,sprK,sprTargetZeroStrain,sprFoldingSequence, ...
        creaseRef,oldCreaseNum,assembleConstant,...
        modelGeometryConstant,modelMechanicalConstant,...
        supportInfo,load);
    
    % Record the diformation field
    if i==1
        UhisThermal(i,:,:)=U;
    else
        UhisThermal(i,:,:)=U+squeeze(UhisThermal(i-1,:,:));
    end
    energyHisThermal(i,:)=squeeze(...
        strainEnergyAssemble(assembleConstant(1),:));    
    tempNode=tempNode+U; 

    % Solve the fold angle for plotting
    node1=squeeze(gravityShape(1,:))'+squeeze(UhisThermal(i,1,:));
    node2=squeeze(gravityShape(2,:))'+squeeze(UhisThermal(i,2,:));
    node3=squeeze(gravityShape(5,:))'+squeeze(UhisThermal(i,5,:));
    node4=squeeze(gravityShape(6,:))'+squeeze(UhisThermal(i,6,:));
    vec1=node2-node1;
    vec2=node4-node3;
    rotation=dot(vec1,vec2)/norm(vec1)/norm(vec2);
    rotation=sign(node4(3)-node1(3))*acos(rotation);
    foldHis(i)=rotation;
    
    % update the thermal matrix for next step
    [thermalMat,thermalNodeNum]=Thermal_AssembleConductMat...
        (gravityShape,newPanel,barConnect,barLength,barType,...
        squeeze(UhisThermal(i,:,:)),modelThermalConstant,...
        modelMechanicalConstant,newCrease2OldCrease,newPanel2OldPanel);
end

assembleNode=tempNode;
Plot_DeformedShapeTemp(viewControl,newNode,assembleNode,newPanel...
    ,squeeze(temperatureHistory(:,thermalStep)),...
    thermalBoundaryPanels,newPanel2OldPanel);
Plot_DeformedHisTemp(viewControl,gravityShape,...
    UhisThermal,newPanel,temperatureHistory);
figure; plot(qhis,foldHis);
figure; plot(foldHis,tempHis);

