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
N=8;
t=0.025;
f=100000; 
connectorBarArea=0.000258;
hingeStiff=0.000001;

gap=0.01;
% Stiffening factor for locking panels to create Yoshimura Folding mode

% Directly Generate the structure using the following function
ori=GenerateModularOrigami(L,N,t,f,gap);
ori.plotBars=1;
ori.showNumber=1;
ori.plotUndeformedShape=0;


% Plot the results for inspection
ori.viewAngle1=145;
ori.viewAngle2=45;
ori.displayRange=[-0.5,1.2,-0.2,2.2,-0.5,1]'; % plotting range
ori.displayRangeRatio=0.1; % plotting range in the negative axis

% plot for inspection
ori.Plot_MeshedOrigami(); % Plot the meshed origami for inspection;

%% Assign Mechanical Properties
% Set up Young's Moduli
ori.panelE=3.18*10^9; 

% ori.creaseE=4*10^9; 
% ori.panelPoisson=0.3;
% ori.creasePoisson=0.3; 


index=[1:552]';
ori.continuingLoading=1;
%% Loading to trigger Yoshimura pattern

nr=ControllerNRLoading;
nr.supp=[2,1,1,1;
         49,0,0,1;
         45,0,0,1;
         92,0,1,1;
         98,1,0,1;
         145,0,0,1;
         141,0,0,1;
         188,0,0,1;
         194,1,0,1;
         241,0,0,1;         
         237,0,0,1;
         284,0,0,1;
         290,1,0,1;
         337,0,0,1;         
         333,0,0,1;
         380,0,0,1;];

loadForce=5*10^-8;
nr.load=cat(2,index,zeros(552,2),-loadForce*ones(552,1));
nr.increStep=50;
nr.tol=1*10^-5;
nr.iterMax=50;
nr.videoOpen=0;


%% Loading to completely fold Yoshimura Pattern

dc=ControllerDCLoading;
dc.supp=[2,1,1,1;
         49,0,0,1;
         45,0,0,1;
         92,0,1,1;
         98,1,0,1;
         145,0,0,1;
         141,0,0,1;
         188,0,0,1;
         194,1,0,1;
         241,0,0,1;         
         237,0,0,1;
         284,0,0,1;
         290,1,0,1;
         337,0,0,1;         
         333,0,0,1;
         380,0,0,1;];

loadForce2=2*10^-5;
dc.load=[45,-loadForce2,0,0;
          92,-loadForce2,0,0;
          141,-loadForce2,-loadForce2,0;
          188,-loadForce2,-loadForce2,0;
          237,-loadForce2,-loadForce2,0;
          284,-loadForce2,-loadForce2,0;
          333,-loadForce2,-loadForce2,0;
          380,-loadForce2,-loadForce2,0;
          369,0,-loadForce2,0;
          374,0,-loadForce2,0;
          345,0,-loadForce2,0;
          350,0,-loadForce2,0;]; 

dc.increStep=45;
dc.tol=5*10^-6;
dc.iterMax=50;
dc.selectedRefDisp=[350,2];
dc.videoOpen=0;

ori.loadingController{1}={"NR",nr};
ori.loadingController{2}={"DC",dc};
ori.Solver_Solve()


%% Reverse the loading for unfolding back to flat state

ori.continuingLoading=1;
ori.loadingController={};

dc_un=ControllerDCLoading;
dc_un.supp=[2,1,1,1;
         49,0,0,1;
         45,0,0,1;
         92,0,1,1;
         98,1,0,1;
         145,0,0,1;
         141,0,0,1;
         188,0,0,1;
         194,1,0,1;
         241,0,0,1;         
         237,0,0,1;
         284,0,0,1;
         290,1,0,1;
         337,0,0,1;         
         333,0,0,1;
         380,0,0,1;];

dc_un.increStep=65;
dc_un.tol=5*10^-6;
dc_un.iterMax=50;
dc_un.selectedRefDisp=[380,1];
dc_un.videoOpen=0;

loadForce2=-8*10^-5;
dc_un.load=[45,-loadForce2,0,0;
          92,-loadForce2,0,0;
          141,-loadForce2,-loadForce2,0;
          188,-loadForce2,-loadForce2,0;
          237,-loadForce2,-loadForce2,0;
          284,-loadForce2,-loadForce2,0;
          333,-loadForce2,-loadForce2,0;
          380,-loadForce2,-loadForce2,0;
          369,0,-loadForce2,0;
          374,0,-loadForce2,0;
          345,0,-loadForce2,0;
          350,0,-loadForce2,0;]; 

ori.loadingController{1}={"DC",dc_un};
ori.Solver_Solve()

nr_un=ControllerNRLoading;
nr_un.supp=[2,1,1,1;
         49,0,0,1;
         45,0,0,1;
         92,0,1,1;
         98,1,0,1;
         145,0,0,1;
         141,0,0,1;
         188,0,0,1;
         194,1,0,1;
         241,0,0,1;
         237,0,0,1;
         284,0,0,1;
         290,1,0,1;
         337,0,0,1;
         333,0,0,1;
         380,0,0,1;];

nr_un.load=cat(2,index,-ori.currentAppliedForce/50);
nr_un.increStep=50;
nr_un.tol=1*10^-5;
nr_un.iterMax=100;
nr_un.videoOpen=0;

ori.loadingController{1}={"NR",nr_un};
ori.Solver_Solve()

%% Folding to form the bridge shape
% First need to adjust the crease stiffness to simulate locking creases

for i=1:N
    ori.sprK(785+119*(i-1)+17)=ori.sprK(785+119*(i-1)+17)/f;
    ori.sprK(785+119*(i-1)+17*3)=ori.sprK(785+119*(i-1)+17*3)/f;
    ori.sprK(785+119*(i-1)+17*5)=ori.sprK(785+119*(i-1)+17*5)/f;
end

for i=1:N
    ori.sprK(785+119*(i-1))=ori.sprK(785+119*(i-1))*f;
    ori.sprK(785+119*(i-1)+17*2)=ori.sprK(785+119*(i-1)+17*2)*f;
    ori.sprK(785+119*(i-1)+17*3)=ori.sprK(785+119*(i-1)+17*3)*f;
    ori.sprK(785+119*(i-1)+17*4)=ori.sprK(785+119*(i-1)+17*4)*f;
    ori.sprK(785+119*(i-1)+17*6)=ori.sprK(785+119*(i-1)+17*6)*f;
end

ori.sprK(1720:2196)=ori.sprK(785+119*(i-1)+17*6)*f;

nr_br=ControllerNRLoading;
nr_br.supp=[2,1,1,1;
         49,0,0,1;
         45,0,0,1;
         92,0,1,1;
         98,1,0,1;
         145,0,0,1;
         141,0,0,1;
         188,0,0,1;
         194,1,0,1;
         241,0,0,1;
         237,0,0,1;
         284,0,0,1;
         290,1,0,1;
         337,0,0,1;
         333,0,0,1;
         380,0,0,1;];

loadForce=10^-9;
nr_br.load=cat(2,index,zeros(552,2),-loadForce*ones(552,1));
nr_br.increStep=60;
nr_br.tol=1*10^-6;
nr_br.iterMax=50;
nr_br.videoOpen=0;

ori.loadingController{1}={"NR",nr_br};
ori.Solver_Solve()

% A two step folding is used to fold to the bridge shape
index2=find(ori.newNode(:,1)>0.95);     
A=size(index2);

nr_br2=ControllerNRLoading;
nr_br2.supp=[2,1,1,1;
         49,0,0,1;
         45,0,0,1;
         92,0,1,1;
         98,1,0,1;
         145,0,0,1;
         141,0,0,1;
         188,0,0,1;
         194,1,0,1;
         241,0,0,1;
         237,0,0,1;
         284,0,0,1;
         290,1,0,1;
         337,0,0,1;
         333,0,0,1;
         380,0,0,1;];

loadForce=2*10^-8;
nr_br2.load=cat(2,index2,-loadForce*ones(A(1),1),zeros(A(1),2));
nr_br2.increStep=37;
nr_br2.tol=1*10^-6;
nr_br2.iterMax=50;
nr_br2.videoOpen=0;

ori.loadingController{1}={"NR",nr_br2};
ori.Solver_Solve()

%% Loading of origami

% Release the strain energy in folding creases
ori.currentSprZeroStrain=ori.Spr_Theta(ori.currentU,ori.sprIJKL,ori.newNode);
ori.currentAppliedForce=zeros(size(ori.newNode));
ori.newNode=ori.newNode+ori.currentU;

% add additional connecting structure to simulate the assembled and
% stiffened structure
inPlateBarNum=96*N;
for i=1:N    
    for j=1:4    
        if (mod(i,2)==1 && mod(j,2)==1) || (mod(i,2)==0 && mod(j,2)==0)
            
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+1,...
                                      48*(i-1)+12*(j-1)+7,...
                                      48*(i-1)+12*(j-1)+3,...
                                      48*(i-1)+12*(j-1)+9,...
                                      48*(i-1)+12*(j-1)+4, ...
                                      48*(i-1)+12*(j-1)+10,...
                                      48*(i-1)+12*(j-1)+6,...
                                      48*(i-1)+12*(j-1)+12,...
                                      connectorBarArea/2,connectorBarArea/5,L,gap,t,hingeStiff)
     
            if j~=4
                ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+8,...
                          48*(i-1)+12*(j-1)+13,...
                          48*(i-1)+12*(j-1)+9,...
                          48*(i-1)+12*(j-1)+14,...
                          48*(i-1)+12*(j-1)+11, ...
                          48*(i-1)+12*(j-1)+16,...
                          48*(i-1)+12*(j-1)+12,...
                          48*(i-1)+12*(j-1)+17,...
                          connectorBarArea/2,connectorBarArea/5,L,gap,t,hingeStiff)
                
            end        
        else
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+3,...
                          48*(i-1)+12*(j-1)+8,...
                          48*(i-1)+12*(j-1)+2,...
                          48*(i-1)+12*(j-1)+7,...
                          48*(i-1)+12*(j-1)+6, ...
                          48*(i-1)+12*(j-1)+11,...
                          48*(i-1)+12*(j-1)+5,...
                          48*(i-1)+12*(j-1)+10,...
                          connectorBarArea/2,connectorBarArea/5,L,gap,t,hingeStiff)

            if j~=4
                
                ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+8,...
                          48*(i-1)+12*(j-1)+13,...
                          48*(i-1)+12*(j-1)+9,...
                          48*(i-1)+12*(j-1)+14,...
                          48*(i-1)+12*(j-1)+11, ...
                          48*(i-1)+12*(j-1)+16,...
                          48*(i-1)+12*(j-1)+12,...
                          48*(i-1)+12*(j-1)+17,...
                          connectorBarArea/2,connectorBarArea/5,L,gap,t,hingeStiff)

            end
        end
    end    
end

%% Here we start adding the bottom connector
btmConnectorNum=inPlateBarNum+49*N;
for i=1:N-1
    for j=1:4
        if (mod(i,2)==1 && mod(j,2)==1) || (mod(i,2)==0 && mod(j,2)==0)
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+6,...
              48*(i-1)+12*(j-1)+54,...
              48*(i-1)+12*(j-1)+5,...
              48*(i-1)+12*(j-1)+52,...
              48*(i-1)+12*(j-1)+3, ...
              48*(i-1)+12*(j-1)+51,...
              48*(i-1)+12*(j-1)+2,...
              48*(i-1)+12*(j-1)+49,...
              connectorBarArea/2,connectorBarArea/5,L,gap,t,hingeStiff)
            
        else
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+10,...
              48*(i-1)+12*(j-1)+58,...
              48*(i-1)+12*(j-1)+12,...
              48*(i-1)+12*(j-1)+59,...
              48*(i-1)+12*(j-1)+7, ...
              48*(i-1)+12*(j-1)+55,...
              48*(i-1)+12*(j-1)+9,...
              48*(i-1)+12*(j-1)+56,...
              connectorBarArea/2,connectorBarArea/5,L,gap,t,hingeStiff)

        end
    end
end

for i=1:N
    ori.sprK(785+119*(i-1))=ori.sprK(785+119*(i-1))/f;
    ori.sprK(785+119*(i-1)+17*2)=ori.sprK(785+119*(i-1)+17*2)/f;
    ori.sprK(785+119*(i-1)+17*3)=ori.sprK(785+119*(i-1)+17*3)/f;
    ori.sprK(785+119*(i-1)+17*4)=ori.sprK(785+119*(i-1)+17*4)/f;
    ori.sprK(785+119*(i-1)+17*6)=ori.sprK(785+119*(i-1)+17*6)/f;
end
ori.sprK(1720:2196)=ori.sprK(785+119*(i-1)+17*6)/f;
ori.sprK=ori.sprK*10^7;

newNodeNum = size(ori.newNode);
newNodeNum = newNodeNum(1);
ori.currentAppliedForce = zeros(newNodeNum,3);   
ori.currentU = zeros(newNodeNum,3);
ori.currentSprZeroStrain = ori.Spr_Theta(ori.currentU,ori.sprIJKL,ori.newNode);


index2=find(ori.newNode(:,1)>0.95);     
A=size(index2);

nr_load=ControllerNRLoading;
nr_load.supp=[37,1,1,1;
              40,1,1,1;
              11,1,1,1;
              8,1,1,1;
              348,1,1,1;
              345,1,1,1;
              374,1,1,1;
              377,1,1,1;];

loadForce=11;

nr_load.load=[175,0,0,-loadForce;
              223,0,0,-loadForce;
              159,0,0,-loadForce;
              207,0,0,-loadForce;
              177,0,0,-loadForce;
              224,0,0,-loadForce;
              158,0,0,-loadForce;
              205,0,0,-loadForce;
              ];
          
          
nr_load.increStep=10;
nr_load.tol=1*10^-5;
nr_load.iterMax=50;
nr_load.videoOpen=0;

ori.loadingController{1}={"NR",nr_load};
ori.Solver_Solve()

ori.Plot_DeformedShape(ori.newNode,ori.newNode+ori.currentU*20);

% Uhis=cat(1,nr.Uhis,dc.Uhis,dc_un.Uhis,nr_un.Uhis,nr_br.Uhis,nr_br2.Uhis );
% ori.Plot_DeformedHis(ori.newNode,Uhis)

ori.plotUndeformedShape=1;
loadIndex=nr_load.load(:,1);
disp=ori.currentU(loadIndex,3);

% Get the displacement-force graph of this loading test
trussIndex=[182,229,183,200];
dispHis=-sum(nr_load.Uhis(:,trussIndex,3),2)*1000/length(ori.currentU(trussIndex,3));
loadHis=(1:nr_load.increStep)*loadForce*8;
figure
plot(dispHis,loadHis);
xlabel('displacement (mm)') 
ylabel('forces (N)') 
grid on
set(gcf,'color','w');
ori.Plot_BarStrain(nr_load)

% These two are the stress values for the most tensile and most compressive
% truss members within the origami bridge
BarStress1=(nr_load.barSxHis(10,460)+nr_load.barSxHis(10,467)+nr_load.barSxHis(10,457))/3
BarStress2=(nr_load.barSxHis(10,473)+nr_load.barSxHis(10,478)+nr_load.barSxHis(10,470))/3

toc