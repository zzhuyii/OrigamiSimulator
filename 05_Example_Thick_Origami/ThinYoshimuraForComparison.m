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

%% Set up the solver
ori=OrigamiSolver;

% build the thin Yoshimura Pattern Manually
ori.panelInnerBarStart=1;
panelBarArea=0.000258;
connectorBarArea=0.000258;

hingeStiff=0.000001;
f=100; %Factor used to stiffen rotational spring to lock them
hingeStiffYoshimura=f*hingeStiff;

% Initialize the size
ori.barArea=[];
ori.barLength=[];
ori.barType=[];
ori.sprK=[];

% Here we will manually add the springs and bars
L=0.25;
N=8;

for i=1:N+1
    ori.newNode(5*(i-1)+1,:)=[L*(i-1),0,0];
    ori.newNode(5*(i-1)+2,:)=[L*(i-1),L,0];
    ori.newNode(5*(i-1)+3,:)=[L*(i-1),2*L,0];
    ori.newNode(5*(i-1)+4,:)=[L*(i-1),3*L,0];
    ori.newNode(5*(i-1)+5,:)=[L*(i-1),4*L,0];
end

ori.AddBar(1,2,panelBarArea,L);
ori.AddBar(2,3,panelBarArea,L);
ori.AddBar(3,4,panelBarArea,L);
ori.AddBar(4,5,panelBarArea,L);

for i=1:N
        ori.AddBar(1+5*(i-1),6+5*(i-1),panelBarArea,L);
        ori.AddBar(5*i,5*(i+1),panelBarArea,L);
    if mod(i,2)==1
        % Diagonal Bars
        ori.AddHinge(1+5*(i-1),2+5*i,2+5*(i-1),1+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        ori.AddHinge(3+5*(i-1),2+5*i,2+5*(i-1),3+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        ori.AddHinge(3+5*(i-1),4+5*i,4+5*(i-1),3+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        ori.AddHinge(5+5*(i-1),4+5*i,4+5*(i-1),5+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        % Horizonal Bars
        ori.AddHinge(1+5*i,2+5*i,1+5*(i-1),1+5*(i+1),...
            panelBarArea,L,hingeStiff,pi);
        ori.AddHinge(2+5*i,3+5*i,3+5*(i-1),3+5*(i+1),...
            panelBarArea,L,hingeStiff,pi);
        ori.AddHinge(3+5*i,4+5*i,3+5*(i-1),3+5*(i+1),...
            panelBarArea,L,hingeStiff,pi);
        ori.AddHinge(4+5*i,5+5*i,5+5*(i-1),5+5*(i+1),...
            panelBarArea,L,hingeStiff,pi);
        % Vertical Bars
        ori.AddHinge(2+5*(i-1),2+5*i,1+5*(i-1),3+5*(i-1),...
            panelBarArea,L,hingeStiffYoshimura,pi);
        ori.AddHinge(3+5*(i-1),3+5*i,2+5*i,4+5*i,...
            panelBarArea,L,hingeStiffYoshimura,pi);
        ori.AddHinge(4+5*(i-1),4+5*i,3+5*(i-1),5+5*(i-1),...
            panelBarArea,L,hingeStiffYoshimura,pi);

    else
        % Diagonal Bars
        ori.AddHinge(2+5*(i-1),1+5*i,1+5*(i-1),2+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        ori.AddHinge(2+5*(i-1),3+5*i,3+5*(i-1),2+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        ori.AddHinge(4+5*(i-1),3+5*i,3+5*(i-1),4+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        ori.AddHinge(4+5*(i-1),5+5*i,5+5*(i-1),4+5*i,...
            panelBarArea,L*sqrt(2),hingeStiff,pi);
        % Horizonal Bars
        if i==N
            ori.AddBar(1+5*i,2+5*i,...
                panelBarArea,L);
            ori.AddBar(2+5*i,3+5*i,...
                panelBarArea,L);
            ori.AddBar(3+5*i,4+5*i,...
                panelBarArea,L);
            ori.AddBar(4+5*i,5+5*i,...
                panelBarArea,L);
        else
            ori.AddHinge(1+5*i,2+5*i,2+5*(i-1),2+5*(i+1),...
                panelBarArea,L,hingeStiff,pi);
            ori.AddHinge(2+5*i,3+5*i,2+5*(i-1),2+5*(i+1),...
                panelBarArea,L,hingeStiff,pi);
            ori.AddHinge(3+5*i,4+5*i,4+5*(i-1),4+5*(i+1),...
                panelBarArea,L,hingeStiff,pi);
            ori.AddHinge(4+5*i,5+5*i,4+5*(i-1),4+5*(i+1),...
                panelBarArea,L,hingeStiff,pi);
        end
        % Vertical Bars
        ori.AddHinge(2+5*(i-1),2+5*i,1+5*i,3+5*i,...
            panelBarArea,L,hingeStiffYoshimura,pi);
        ori.AddHinge(3+5*(i-1),3+5*i,2+5*(i-1),4+5*(i-1),...
            panelBarArea,L,hingeStiffYoshimura,pi);
        ori.AddHinge(4+5*(i-1),4+5*i,3+5*i,5+5*i,...
            panelBarArea,L,hingeStiffYoshimura,pi);
    end
end

%% Set up terms for initialization
newNodeNum = size(ori.newNode);
newNodeNum = newNodeNum(1);
ori.currentAppliedForce = zeros(newNodeNum,3);   
ori.currentU = zeros(newNodeNum,3);
ori.currentSprZeroStrain = pi*logical(ori.sprK);



% Plot the results for inspection
ori.viewAngle1=45;
ori.viewAngle2=45;
ori.displayRange=[-0.5,3,-0.5,3,-1,1]'; % plotting range
ori.plotBars=1;
ori.showNumber=1;

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

LoadStep=20;
eigVal=zeros(LoadStep,15);

nr=ControllerNRLoading;

nr.supp=[49,1,1,1;
         56,0,0,1;
         53,0,1,1;
         60,0,0,1];

loadForce=0.0005;

nr.load=[63,0,0,-loadForce;
         64,0,0,loadForce;
         61,0,0,-loadForce;
         62,0,0,loadForce];

% Mode 2
nr.load=[63,0,0,loadForce;
         64,0,0,loadForce;
         61,0,0,loadForce;
         62,0,0,loadForce];

% Mode 3
LoadStep=11;
loadForce=0.001;
nr.supp=[1,1,1,1;
         39,0,0,1;
         31,0,1,1;
         26,0,0,1];
     
nr.load=[18,0,0,loadForce;
         23,0,0,loadForce;
         16,0,0,loadForce;
         11,0,0,loadForce;
         24,0,0,loadForce;
         47,0,0,loadForce;];

nr.increStep=1;
nr.tol=1*10^-7;
nr.iterMax=100;
nr.videoOpen=0;
nr.plotOpen=0;



for i=1:LoadStep
    
    % Find the eigen shape
    StiffMat=ori.Solver_CalcK();

    [Vec,Val]=eigs(StiffMat,15,'smallestabs');
    Vec=real(Vec);
    Val=real(Val);
    Val=diag(Val);
    
    eigVal(i,:)=abs(Val');
    
    ori.loadingController{1}={"NR",nr};
    ori.Solver_Solve()

end

[Vec,Val]=eigs(StiffMat,15,'smallestabs');
Vec=real(Vec);
Val=real(Val);
Val=diag(Val);

U1=Vec(:,1);
U1=reshape(-U1,3,64)';
ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
    ori.newNode+U1+ori.currentU)

%% Find the mode shape at flat

% StiffMat=ori.Solver_CalcK();
% [Vec,Val]=eigs(StiffMat,15,'smallestabs');
% Vec=real(Vec);
% Val=real(Val);
% Val=diag(Val);
% 
% 
% U7=Vec(:,7);
% U8=Vec(:,8);
% U9=Vec(:,9);
% 
% U7=reshape(-U7,3,64)';
% U8=reshape(-U8,3,64)';
% U9=reshape(-U9,3,64)';
% 
% 
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+U7+ori.currentU)
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+U8+ori.currentU)
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+U9+ori.currentU)
% 
% 
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+(0.8*U7+U8+0.5*U9)+ori.currentU)
% 
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+(U7-U8+0.5*U9)+ori.currentU)
% 
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+(U7-1*U9)+ori.currentU)
% 
% ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
%     ori.newNode+(U7-0.45*U9)+ori.currentU)


toc


