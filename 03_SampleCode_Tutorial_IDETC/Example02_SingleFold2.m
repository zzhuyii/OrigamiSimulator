%%%%% Sequentially Working Origami Multi-Physics Simulator (SWOMPS)  %%%%%%
%
% Authors: Yi Zhu, and Evgueni T. Filipov
%
% Discription: This code package implement a bar and hinge model based 
% simulator for active origami structures with multi-physics based 
% actuation mechanisms. The code package can capture both the mechanical 
% behavior and the heat transfer aspect. The implementation is versatile
% and has the following features:
%
% (1) Provides 5 different loading solvers of active origami. They are: 
%     Newton-Raphson method, displacement controlled method, modified 
%     generazlied displacement controlled method, self-stress folding, and
%     thermal folding method.
% (2) Allows users to create arbitrary number and sequence of the five
%     loading methods. Users can stop the solver at specified increments
%     and switch between different solvers or edit origami systems during 
%     within the increment easily.
% (3) Simulate electro-thermo-mechanically coupled actuation of origami.
% (4) Simulate inter-panel contact of origami systems.
% (5) Simulate the compliant creases explicitly with novel bar and hinge
%     model meshing schemes.
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% origami simulator, the origami contact, compliant crease, electro-thermal
% model presented in this package. 
%
% Reference:
% [1] Y. Zhu, E. T. Filipov (2021). 'Sequentially Working Origami Multi-
%     Physics Simulator (SWOMPS): A Versatile Implementation' (submitted)
% [2] Y. Zhu, E. T. Filipov (2021). 'Rapid Multi-Physic Simulation for 
%     Electro-Thermal Origami Robotic Systems' (submitted)
% [3] Y. Zhu, E. T. Filipov (2020). 'A Bar and Hinge Model for Simulating 
%     Bistability in Origami Structures with Compliant Creases' Journal of 
%     Mechanisms and Robotics, 021110-1. 
% [4] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' Proc. R. Soc. A, 475: 
%     20190366.       
% [5] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [6] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [7] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' Proc. R. Soc. A 473: 
%     20170348. 
% [8] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%% Sequentially Working Origami Multi-Physics Simulator (SWOMPS)  %%%%%%
close all; clear all;

%% Initialize the solver
ori=OrigamiSolver;

%% Input data for Defining the Geometry

% input that will be recorded
w=(rand*300+80)*(10^-6);
l=1000*(10^-6);
t1=(rand*0.2+0.1)*(10^-6);
t2=(rand*0.4+0.6)*(10^-6); 
panelThick=(rand*15+5)*(10^-6);

% other input we will fix
underCut=100*(10^-6);
rho=1200;

plotFlag=1;

ori.node0=[0 0 0;
           l 0 0;
           0 l+0.5*w 0;
           l l+0.5*w 0;
           0 l+l+w 0;
           l l+l+w 0;        
           -l -l -underCut;
           2*l -l -underCut;
           -l 3*l -underCut;
           2*l 3*l -underCut;];

ori.panel0{1}=[1 2 4 3];
ori.panel0{2}=[3 4 6 5];
ori.panel0{3}=[7 8 10 9];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=5*10^(-3); % plotting range
ori.displayRangeRatio=0.5; % plotting range in the negative axis

% Plot the unmeshed origami for inspection;
% ori.Plot_UnmeshedOrigami(); 

%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthMat=zeros(ori.oldCreaseNum,1);
ori.creaseWidthMat(4)=w;

% Compute the meshed geometry
ori.Mesh_CompliantCreaseGeometry()
% Plot the meshed origami for inspection;
% ori.Plot_MeshedOrigami(); 


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickMat=[panelThick;panelThick;panelThick]; 
ori.panelW=w;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=100;

ori.creaseThickMat=zeros(ori.oldCreaseNum,1);
ori.creaseThickMat(4)=(t1+t2);

%% setup panel contact information

ori.contactOpen=1;
ori.ke=0.0001;
ori.d0edge=40*(10^(-6));
ori.d0center=40*(10^(-6));


%% Assign Thermal Properties

ori.panelThermalConductMat = [0.3;0.3;0.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=(l)*1.5;

%% Define Density of Origami
ori.densityCrease=rho;
ori.densityPanel=rho;


%% Setup the loading controller

% applying the gravity loading
nr=ControllerNRLoading;

nr.increStep=80;
nr.tol=2*10^-6;
nr.iterMax=50;

nr.supp=[1,1,1,1;
      2,1,1,1;
      3,1,1,1;
      4,1,1,1;   
      9,1,1,1;
      10,1,1,1;
      11,1,1,1;
      12,1,1,1;];  
  
loadMag=l*l*panelThick*rho/6*9.8;
loadMag=loadMag/nr.increStep;

nr.load=[5,0,0,-loadMag;
      6,0,0,-loadMag;
      7,0,0,-loadMag;
      8,0,0,-loadMag;
      17,0,0,-2*loadMag;];
nr.videoOpen=0;

ori.loadingController{1}={"NR",nr};

%% Solving the model
ori.Solver_Solve();

% record the deformation before thermal loading
gravityU=ori.currentU;


%% Proceed to the thermal loading
% we will continue the loading process
ori.continuingLoading=1;

Qstep=0.2*10^-3;
maxLoadingStep=500;

newNodeNum=size(ori.newNode);
newNodeNum=newNodeNum(1);

tempHisAssemble=zeros(newNodeNum,maxLoadingStep);
UhisAssemble=zeros(maxLoadingStep,newNodeNum,3);

finalQBase=0;
finalTBase=0;
finalStepBase=0;

%% first we load the base creases to assemble

for i=1:maxLoadingStep
    thermal=ControllerThermalLoading;
    thermal.thermalStep=1;
    thermal.tol=5*10^-7; 

    thermal.supp=[1,1,1,1;
                  2,1,1,1;
                  3,1,1,1;
                  4,1,1,1;   
                  9,1,1,1;
                  10,1,1,1;
                  11,1,1,1;
                  12,1,1,1;];  

    thermal.thermalBoundaryPanelMat=[3];
    thermal.roomTempNode=[];

    thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
    thermal.deltaAlpha(4)=(52-14)*10^(-6);
    
    thermal.Emat1=0.5*79*10^9; 
    thermal.Emat2=2*10^9;
    thermal.tmat1=t1;
    thermal.tmat2=t2;
    thermal.videoOpen=0; % close the animation
    thermal.plotOpen=0; % close the plot

    % the target loading of crease heating
    thermal.targetCreaseHeating=[4,Qstep];
    ori.loadingController{1}={"ThermalLoading",thermal};
    ori.Solver_Solve();
    % we perform the countinuing loading step by step and check the
    % rotation angle, the loading stops after reaching 90 degree;
    
    UhisAssemble(i,:,:)=thermal.Uhis(1,:,:);
    tempHisAssemble(:,i)=thermal.temperatureHis(:,1);
    
    % Solve the fold angle for stop loading
    node1=squeeze(ori.newNode(1,:))'+squeeze(UhisAssemble(i,1,:));
    node2=squeeze(ori.newNode(4,:))'+squeeze(UhisAssemble(i,4,:));
    node3=squeeze(ori.newNode(5,:))'+squeeze(UhisAssemble(i,5,:));
    node4=squeeze(ori.newNode(8,:))'+squeeze(UhisAssemble(i,8,:));
    vec1=node2-node1;
    vec2=node4-node3;
    rotation=dot(vec1,vec2)/norm(vec1)/norm(vec2);
    rotation=sign(node4(3)-node1(3))*acos(rotation);
    % Solve the average temperature of crease
    Tave=1/3*(tempHisAssemble(13,i)+...
        tempHisAssemble(14,i)+...
        tempHisAssemble(15,i));
    
    if rotation>pi/2
        finalQBase=Qstep*i;
        finalTBase=Tave;  
        finalStepBase=i;
        break
    end   
end

if plotFlag==1
    % plot the loaded structure after the assemble step
    ori.Plot_DeformedShapeTemp(thermal,...
        gravityU+ori.newNode, ori.currentU+ori.newNode, ...
        squeeze(tempHisAssemble(:,finalStepBase)));
end


%% Analyze the current freqency

frequency=ControllerFrequencyAnalysis;
frequency.supp=[1,1,1,1;
                  2,1,1,1;
                  3,1,1,1;
                  4,1,1,1;   
                  9,1,1,1;
                  10,1,1,1;
                  11,1,1,1;
                  12,1,1,1;]; 
          
[frequencySquared,Umode]=ori.Dynamic_FrequencyAnalysis(frequency);
frequencySquared=diag(frequencySquared);
[freq,index]=min(frequencySquared);
freqFinal=sqrt(freq);
UmodeBase=Umode(:,index);
UmodeBase=reshape(UmodeBase,3,18)'/10000;
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,ori.newNode+ori.currentU+UmodeBase);



fprintf('Simulation Finished\n');
fprintf('Final input energy for assemble: %e W;\n',finalQBase);
fprintf('Final crease temperature for assemble: %e degree C;\n',finalTBase);


outputVector=[w t1 t2 panelThick finalQBase finalTBase freqFinal];

