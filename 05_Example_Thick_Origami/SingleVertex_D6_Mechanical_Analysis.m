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

%% Initialize the solver
clear;clc;close all;
tic

%% Define the Geometry of origami
L=0.25;
t=0.1;
gap=0.1;

% Directly Generate the structure using the following function
ori=GenerateThickYoshimuraUnitCell(L,t,gap);
ori.plotBars=1;
ori.showNumber=1;
ori.plotUndeformedShape=1;
ori.compliantCreaseOpen=0;

% Plot the results for inspection
ori.viewAngle1=245;
ori.viewAngle2=45;
ori.displayRange=[-0.2,1,-0.2,1,-0.3,0.7]'; % plotting range
ori.displayRangeRatio=0.5; % plotting range in the negative axis

% plot for inspection
ori.plotBars=1;
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;

%% Assign Mechanical Properties
% Set up Young's Moduli
ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 


%% Solve for the eigen mode at flat
StiffMat=ori.Solver_CalcK();

[Vec,Val]=eigs(StiffMat,15,'smallestabs');
Vec=real(Vec);
Val=real(Val);
Val=diag(Val);

Umode7=Vec(:,7);
Umode8=Vec(:,8);
Umode9=Vec(:,9);

ori.Plot_DeformedShape(ori.newNode,ori.newNode+0.5*[reshape(Umode7,[3,168/3])]');
ori.Plot_DeformedShape(ori.newNode,ori.newNode+0.5*[reshape(Umode8,[3,168/3])]');
ori.Plot_DeformedShape(ori.newNode,ori.newNode+0.5*[reshape(Umode9,[3,168/3])]');


%% Yoshimura Folding Target
foldTarget=readtable('FoldHis.xlsx');
foldTarget=table2array(foldTarget);
%foldTarget=foldTarget(2:16,:);


%% Solve how the eigen value changes when folding with Yoshimura Mode

ori.continuingLoading=1;
ori.Solver_Solve()
% Skip the auto initialization step

selfFold=ControllerSelfFolding();

selfFold.supp=[4,1,1,1;
               5,1,1,1;
               7,1,1,1;
               8,1,1,1;];

selfFold.increStep=5;
selfFold.tol=1*10^(-3);
selfFold.iterMax=500;
selfFold.targetRotZeroStrain=pi*ones(202,1);

LoadStep=15;
selfFold.videoOpen=0;
ori.plotUndeformedShape=0;
eigVal=zeros(15,LoadStep);

index=[117,134,151,168,185,202];
RotationHis=zeros(LoadStep+1,6);
RotationHis(1,:)=pi*ones(1,6);

StiffMat=ori.Solver_CalcK();

[Vec,Val]=eigs(StiffMat,15,'smallestabs');
Vec=real(Vec);
Val=real(Val);
Val=diag(Val);

eigVal(:,1)=abs(Val);


%% Incrementally folding the unit cell 
% We will intcrementally folding the unit cell while tracking the folding
% angle and the eigen vlaues of the stiffness matrix

for i=1:LoadStep
    
    selfFold.targetRotZeroStrain(134)=pi-foldTarget(i,1);
    selfFold.targetRotZeroStrain(117)=pi-foldTarget(i,2);

    selfFold.targetRotZeroStrain(151)=pi-foldTarget(i,3);
    selfFold.targetRotZeroStrain(168)=pi-foldTarget(i,4);   

    selfFold.targetRotZeroStrain(185)=pi+foldTarget(i,5);
    selfFold.targetRotZeroStrain(202)=pi+foldTarget(i,6);

    % Find the eigen shape
    selfFold.plotOpen=0;
    if i==LoadStep
        selfFold.plotOpen=1;
    end

    ori.loadingController{1}={"SelfFold",selfFold};
    ori.Solver_Solve();    

    RotationHis(i+1,:)=selfFold.sprRotHis(5,index);

    StiffMat=ori.Solver_CalcK();

    [Vec,Val]=eigs(StiffMat,15,'smallestabs');
    Vec=real(Vec);
    Val=real(Val);
    Val=diag(Val);

    eigVal(:,i+1)=abs(Val);
end

eigVal=eigVal';

%% Plot the folding curve
FoldHis=-(RotationHis-pi);
FoldHis(:,5:6)=-FoldHis(:,5:6);

figure
hold on
plot3(FoldHis(:,1),FoldHis(:,3),FoldHis(:,5))

set(gcf, 'color', 'white');
set(gca,'DataAspectRatio',[1 1 1])
view(60,15); 
axis([0 4 0 4 0 4]);
hAxis = gca;
hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
hAxis.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
hAxis.YRuler.SecondCrossoverValue  = 0; % Y crossover with Z axis
hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis;
hAxis.FontSize = 14;
xticks([0 pi/2 pi]);
xticklabels({'','',''})
yticks([pi/2 pi]);
yticklabels({'',''})
zticks([pi/2 pi]);
zticklabels({'',''})