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
N=2;
t=0.025;
f=100000; 
gap=0.01;
% Stiffening factor for locking selected creases

% Directly Generate the structure using the following function
ori=GenerateModularOrigami(L,N,t,f,gap);
ori.plotBars=1;
ori.showNumber=0;
ori.plotUndeformedShape=0;


% Plot the results for inspection
ori.viewAngle1=145;
ori.viewAngle2=45;
ori.displayRange=[-0.5,1.5,-0.5,1,-0.5,1]'; % plotting range

% plot for inspection
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;

%% Assign Mechanical Properties
% Set up Young's Moduli
ori.panelE=2*10^9; 
ori.creaseE=2*10^9; 
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

%% Self fold to the configuration
% Skip the auto initialization step
ori.continuingLoading=1;
nr=ControllerNRLoading;
     
nr.supp=[10,1,1,1;
         90,0,0,1;
         58,1,0,0;];      

loadForce=0.000000001;

index=[1:96]';
nr.load=cat(2,index,zeros(96,2),-loadForce*ones(96,1));

nr.increStep=200;
nr.tol=5*10^-6;
nr.iterMax=50;
nr.videoOpen=0;


ori.loadingController{1}={"NR",nr};
ori.Solver_Solve()

for i=1:nr.increStep
    dispHis(i)=mean(nr.Uhis(i,1:48*2,3)-nr.Uhis(i,10,3));
    % Change the second index between 2 and 10 for different support
    % configuration of the system
end
dispHis=dispHis';
% This vector stores the motion of the center of the gravity relative to
% supporting nodes

Fhis=(1:nr.increStep)*(48*2*loadForce);
Fhis=Fhis';
% This vector stores the total gravity applied


% Calculate the internal strain energy from the rotational springs
rotK=ori.sprK';
rotK(rotK == 0)=[];

% lockIndex=[2,4,6,9,11,13]; % Index of locked Creases
% rotK(lockIndex)=rotK(lockIndex)*f;

rotHis=nr.sprRotHis;
rotHis(:,all(rotHis == 0))=[]; % removes column if the entire column is zero
rotE=zeros(nr.increStep,1);

for i=1:nr.increStep
    rotE(i)=sum(rotK.*(rotHis(i,:)-pi).*(rotHis(i,:)-pi))/2;
end
toc