%%%%% Sequentially Working Origami Multi-Physics Simulator (SWOMPS)  %%%%%%
%
% Authors: Yi Zhu, and Evegueni T. Filipov
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

function [finalQBase,finalQFunc,finalTBase,finalTFunc]=GripperFunc(w,la,l1,l2,t1,t2,plotFlag)


%% Initialize the solver
ori=OrigamiSolver;

%% Input data for Defining the Geometry

% w=300*(10^-6);
% la=1000*(10^-6);
% l1=1000*(10^-6);
% l2=1000*(10^-6);
% t1=0.2*(10^-6);
% t2=0.8*(10^-6); 
% plotFlag=1;

range=1000*(10^-6);
underCut=150*(10^-6);
rho=1200;
tpanel=5*(10^-6);

ori.node0=[0 0 0;
      l1 0 0;
      l1+range 0 0;
      2*l1+range 0 0;
      0 l2 0;
      l1 l2 0;
      l1+range l2 0;
      2*l1+range l2 0;
      0.5*l1 la+l2 0;
      l1 la+l2 0;
      l1+range la+l2 0;
      1.5*l1+range la+l2 0;
      -l1 -l1 -underCut;
      3*l1+range -0.5*l1 -underCut;
      -0.5*l1 la+l2+l1 -underCut;
      3*l1+range la+l2+l1 -underCut;];

ori.panel0{1}=[1 2 6 5];
ori.panel0{2}=[2 3 7 6];
ori.panel0{3}=[3 4 8 7];
ori.panel0{4}=[5 6 10 9];
ori.panel0{5}=[7 8 12 11];
ori.panel0{6}=[13 14 16 15];

% Analyze the original pattern before proceeding to the next step
ori.Mesh_AnalyzeOriginalPattern();


%% Meshing of the origami model

% Define the crease width 
ori.creaseWidthVec=zeros(ori.oldCreaseNum,1);
ori.creaseWidthVec(3)=w;
ori.creaseWidthVec(4)=w;
ori.creaseWidthVec(6)=w;
ori.creaseWidthVec(10)=w;

% Compute the meshed geometry
ori.Mesh_Mesh()

% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=5*10^(-3); % plotting range
ori.displayRangeRatio=0.2; % plotting range in the negative axis

ori.newNode(14,:)=ori.newNode(14,:)+[-w/2,0,0];
ori.newNode(15,:)=ori.newNode(15,:)+[-w/2,0,0];
ori.newNode(20,:)=ori.newNode(20,:)+[w/2,0,0];
ori.newNode(17,:)=ori.newNode(17,:)+[w/2,0,0];
ori.newNode(7,:)=ori.newNode(7,:)+[0,-w/2,0];
ori.newNode(8,:)=ori.newNode(8,:)+[0,-w/2,0];
ori.newNode(32,:)=ori.newNode(32,:)+[0,-w/4,0];
ori.newNode(26,:)=ori.newNode(26,:)+[0,-w/4,0];
ori.newNode(38,:)=ori.newNode(38,:)+[0,-w/4,0];
ori.newNode(28,:)=ori.newNode(28,:)+[-w/4,0,0];
ori.newNode(35,:)=ori.newNode(35,:)+[w/4,0,0];
ori.newNode(40,:)=ori.newNode(40,:)+[-w/8,0,0];
ori.newNode(41,:)=ori.newNode(41,:)+[w/8,0,0];

%ori.Plot_UnmeshedOrigami(); % Plot the unmeshed origami for inspection;
%ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;


%% Assign Mechanical Properties

ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 
ori.panelThickVec=[tpanel;tpanel;tpanel;tpanel;tpanel;500*10^(-6)]; 
ori.panelW=w;

% set up the diagonal rate to be large to suppress crease torsion
ori.diagonalRate=100;

ori.creaseThickVec=zeros(ori.oldCreaseNum,1);
ori.creaseThickVec(3)=(t1+t2);
ori.creaseThickVec(4)=(t1+t2);
ori.creaseThickVec(6)=(t1+t2);
ori.creaseThickVec(10)=(t1+t2);

%% setup panel contact information

ori.contactOpen=1;
ori.ke=0.0001;
ori.d0edge=40*(10^(-6));
ori.d0center=40*(10^(-6));


%% Assign Thermal Properties

ori.panelThermalConductVec = [0.3;0.3;0.3;0.3;0.3;1.3]; 
ori.creaseThermalConduct=0.3;
ori.envThermalConduct=0.026;

% thickness of the submerged environment at RT
ori.t2RT=(l1+l2)*1.5/2;



%% Setup the loading controller

% applying the gravity loading
nr=ControllerNRLoading;

nr.increStep=80;
nr.tol=2*10^-6;
nr.iterMax=50;

nr.supp=[5,1,1,1;
      6,1,1,1;
      7,1,1,1;
      8,1,1,1;
      38,1,1,1;      
      21,1,1,1;
      22,1,1,1;
      23,1,1,1;
      24,1,1,1;
      42,1,1,1;];  
  
loadMagBase=l1*l2*tpanel*rho/6*9.8;
loadMagAct=3/4*l1*la*tpanel*rho/6*9.8;

loadMagBase=loadMagBase/nr.increStep;
loadMagAct=loadMagAct/nr.increStep;

nr.load=[1,0,0,-loadMagBase;
      2,0,0,-loadMagBase;
      3,0,0,-loadMagBase;
      4,0,0,-loadMagBase;
      37,0,0,-2*loadMagBase;
      9,0,0,-loadMagBase;
      10,0,0,-loadMagBase;
      11,0,0,-loadMagBase;
      12,0,0,-loadMagBase;
      39,0,0,-2*loadMagBase;
      13,0,0,-loadMagAct;
      14,0,0,-loadMagAct;
      15,0,0,-loadMagAct;
      16,0,0,-loadMagAct;
      40,0,0,-2*loadMagAct;
      17,0,0,-loadMagAct;
      18,0,0,-loadMagAct;
      19,0,0,-loadMagAct;
      20,0,0,-loadMagAct;
      41,0,0,-2*loadMagAct;];
nr.videoOpen=0;
nr.plotOpen=0;

ori.loadingController{1}={"NR",nr};

%% Solving the model
ori.Solver_Solve();

% record the deformation before thermal loading
gravityU=ori.currentU;


%% Proceed to the thermal loading
% we will continue the loading process
ori.continuingLoading=1;

Qstep=0.1*10^-3;
maxLoadingStep=400;

newNodeNum=size(ori.newNode);
newNodeNum=newNodeNum(1);

tempHisAssemble=zeros(newNodeNum,maxLoadingStep);
UhisAssemble=zeros(maxLoadingStep,newNodeNum,3);

finalQBase=0;
finalTBase=0;
finalStepBase=0;

%% first we load the base creases to assemble

for i=1:maxLoadingStep
    thermal=ControllerElectroThermalFolding;
    thermal.thermalStep=1;
    thermal.tol=5*10^-7; 

    thermal.supp=[5,1,1,1;
              6,1,1,1;
              7,1,1,1;
              8,1,1,1;
              38,1,1,1;      
              21,1,1,1;
              22,1,1,1;
              23,1,1,1;
              24,1,1,1;
              42,1,1,1;];

    thermal.thermalBoundaryPanelVec=[6];
    thermal.roomTempNode=[38];

    thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
    thermal.deltaAlpha(3)=(52-14)*10^(-6);
    thermal.deltaAlpha(4)=(52-14)*10^(-6);
    thermal.deltaAlpha(6)=(52-14)*10^(-6);
    thermal.deltaAlpha(10)=(52-14)*10^(-6);
    
    thermal.Emat1=0.5*79*10^9; 
    thermal.Emat2=2*10^9;
    thermal.tmat1=t1;
    thermal.tmat2=t2;
    thermal.videoOpen=0; % close the animation
    thermal.plotOpen=0; % close the plot

    % the target loading of crease heating
    thermal.targetCreaseHeating=[3,Qstep;6,Qstep];
    ori.loadingController{1}={"ThermalLoading",thermal};
    ori.Solver_Solve();
    % we perform the countinuing loading step by step and check the
    % rotation angle, the loading stops after reaching 90 degree;
    
    UhisAssemble(i,:,:)=thermal.Uhis(1,:,:);
    tempHisAssemble(:,i)=thermal.temperatureHis(:,1);
    
    % Solve the fold angle for stop loading
    node1=squeeze(ori.newNode(1,:))'+squeeze(UhisAssemble(i,1,:));
    node2=squeeze(ori.newNode(2,:))'+squeeze(UhisAssemble(i,2,:));
    node3=squeeze(ori.newNode(5,:))'+squeeze(UhisAssemble(i,5,:));
    node4=squeeze(ori.newNode(6,:))'+squeeze(UhisAssemble(i,6,:));
    vec1=node2-node1;
    vec2=node4-node3;
    rotation=dot(vec1,vec2)/norm(vec1)/norm(vec2);
    rotation=sign(node1(3)-node4(3))*acos(rotation);
    % Solve the average temperature of crease
    Tave1=1/3*(tempHisAssemble(25,i)+...
        tempHisAssemble(26,i)+...
        tempHisAssemble(27,i));

    Tave2=1/3*(tempHisAssemble(31,i)+...
        tempHisAssemble(32,i)+...
        tempHisAssemble(33,i));
    
    if rotation>pi/2
        finalQBase=Qstep*i;
        finalTBase=max(Tave1,Tave2);  
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

tempHisGrip=zeros(newNodeNum,maxLoadingStep);
UhisGrip=zeros(maxLoadingStep,newNodeNum,3);

finalQFunc=0;
finalTFunc=0;
finalStepFunc=0;

%% Next we load the two gripper panel to close the gripping arms
for i=1:maxLoadingStep
    thermal=ControllerElectroThermalFolding;
    thermal.thermalStep=1;
    thermal.tol=5*10^-7; 

    thermal.supp=[5,1,1,1;
              6,1,1,1;
              7,1,1,1;
              8,1,1,1;
              38,1,1,1;      
              21,1,1,1;
              22,1,1,1;
              23,1,1,1;
              24,1,1,1;
              42,1,1,1;];

    thermal.thermalBoundaryPanelVec=[6];
    thermal.roomTempNode=[38];

    thermal.deltaAlpha=zeros(ori.oldCreaseNum,1);
    thermal.deltaAlpha(3)=(52-14)*10^(-6);
    thermal.deltaAlpha(4)=(52-14)*10^(-6);
    thermal.deltaAlpha(6)=(52-14)*10^(-6);
    thermal.deltaAlpha(10)=(52-14)*10^(-6);
    
    thermal.Emat1=0.5*79*10^9; 
    thermal.Emat2=2*10^9;
    thermal.tmat1=t1;
    thermal.tmat2=t2;
    thermal.videoOpen=0; % close the animation
    thermal.plotOpen=0; % close the plot

    % the target loading of crease heating
    thermal.targetCreaseHeating=[4,Qstep;10,Qstep];
    ori.loadingController{1}={"ThermalLoading",thermal};
    ori.Solver_Solve();
    % we perform the countinuing loading step by step and check the
    % rotation angle, the loading stops after reaching 90 degree;
    
    UhisGrip(i,:,:)=thermal.Uhis(1,:,:);
    tempHisGrip(:,i)=thermal.temperatureHis(:,1);
    
    % Solve the distance for stoping when gripper closed
    node1=squeeze(ori.newNode(16,:))'+squeeze(UhisGrip(i,16,:));
    node2=squeeze(ori.newNode(19,:))'+squeeze(UhisGrip(i,19,:));
    disp=norm(node1-node2);
    
    % Solve the average temperature of crease
    Tave1=1/3*(tempHisGrip(34,i)+...
        tempHisGrip(35,i)+...
        tempHisGrip(36,i));

    Tave2=1/3*(tempHisGrip(28,i)+...
        tempHisGrip(29,i)+...
        tempHisGrip(30,i));
    
    if disp<100*(10^-6)
        finalQFunc=Qstep*i;
        finalTFunc=max(Tave1,Tave2);  
        finalStepFunc=i;
        break
    end   
end

if plotFlag==1
    % plot the loaded structure after the assemble step
    ori.Plot_DeformedShapeTemp(thermal,...
        gravityU+ori.newNode, ori.currentU+ori.newNode, ...
        squeeze(tempHisGrip(:,finalStepFunc)));
end

fprintf('Simulation Finished\n');
fprintf('Final input energy for assemble: %e W;\n',finalQBase);
fprintf('Final input energy for gripping: %e W;\n',finalQFunc);
fprintf('Final crease temperature for assemble: %e degree C;\n',finalTBase);
fprintf('Final crease temperature for gripping: %e degree C;\n',finalTFunc);

end