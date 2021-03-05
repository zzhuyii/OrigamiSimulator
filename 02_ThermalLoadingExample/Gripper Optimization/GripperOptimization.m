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

%% This code applies the coodinate descend to optimize the power 
%  consumption of the gripper
tic

OptStep=15;
OptStepRan=0;
tg=0.2*10^-6;
ts=0.8*10^-6;

% starting value of parameters
Wstart=200*10^-6;
LaStart=700*10^-6;
L1Start=1000*10^-6;
L2Start=1000*10^-6;

% step length for each parameters
dW=10*10^-6;
dLa=50*10^-6;
dL1=50*10^-6;
dL2=50*10^-6;

% Process
Qhis=zeros(OptStep,1);
Whis=zeros(OptStep,1);
This=zeros(OptStep,1);
Lahis=zeros(OptStep,1);
L1his=zeros(OptStep,1);
L2his=zeros(OptStep,1);

% optimization range
maxW=300*10^-6;
minW=100*10^-6;
maxLa=1500*10^-6;
minLa=500*10^-6;
maxL1=1500*10^-6;
minL1=500*10^-6;
maxL2=1500*10^-6;
minL2=500*10^-6;

TgSU8=210;
Tmax=1.2*TgSU8;

W=Wstart;
La=LaStart;
L1=L1Start;
L2=L2Start;

[FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1,L2,tg,ts,0);
tempPower=FinalQBase+FinalQFunc;

Qstart=tempPower;
Tstart=max(FinalTBase,FinalTFunc);

for i=1:OptStep
    move=0;
    
    [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W-dW,La,L1,L2,tg,ts,0);
    tempPowerMinus=FinalQBase+FinalQFunc;
    if tempPowerMinus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && W+dW>=minW
        W=W-dW;
        tempPower=tempPowerMinus;
        move=1;
    else
        [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W+dW,La,L1,L2,tg,ts,0);
        tempPowerPlus=FinalQBase+FinalQFunc;
        if tempPowerPlus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && W+dW<=maxW
            W=W+dW;
            tempPower=tempPowerPlus;
            move=1;
        end
    end
        
    [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La-dLa,L1,L2,tg,ts,0);
    tempPowerMinus=FinalQBase+FinalQFunc;
    if tempPowerMinus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && La+dLa>=minLa
        La=La-dLa;
        tempPower=tempPowerMinus;
        move=1;
    else
        [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La+dLa,L1,L2,tg,ts,0);
        tempPowerPlus=FinalQBase+FinalQFunc;
        if tempPowerPlus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && La+dLa<=maxLa
            La=La+dLa;
            tempPower=tempPowerPlus;
            move=1;
        end
    end
    
    [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1-dL1,L2,tg,ts,0);
    tempPowerMinus=FinalQBase+FinalQFunc;
    if tempPowerMinus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && L1-dL1>=minL1
        L1=L1-dL1;
        tempPower=tempPowerMinus;
        move=1;
    else
        [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1+dL1,L2,tg,ts,0);
        tempPowerPlus=FinalQBase+FinalQFunc;
        if tempPowerPlus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && L1+dL1<=maxL1
            L1=L1+dL1;
            tempPower=tempPowerPlus;
            move=1;
        end
    end
    
    [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1,L2-dL2,tg,ts,0);
    tempPowerMinus=FinalQBase+FinalQFunc;
    if tempPowerMinus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && L2-dL2>=minL2
        L2=L2-dL2;
        tempPower=tempPowerMinus;
        move=1;
    else
        [FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1,L2+dL2,tg,ts,0);
        tempPowerPlus=FinalQBase+FinalQFunc;
        if tempPowerPlus<tempPower && FinalTBase<Tmax && FinalTFunc<Tmax && L2+dL2<=maxL2 
            L2=L2+dL2;
            tempPower=tempPowerPlus;
            move=1;
        end
    end
    
    
    Qhis(i)=tempPower;
    Whis(i)=W;
    This(i)=max(FinalTBase,FinalTFunc);
    Lahis(i)=La;
    L1his(i)=L1;
    L2his(i)=L2;
    
    
    OptStepRan=i;
    if move==0
       OptStepRan=1;
       fprintf('Optimization Stopped after %d cycles\n',i);
       i=OptStep+1; 
    end
end
t=toc;

[StartlQBase,StartQFunc,StartTBase,StartFunc]=GripperFunc(Wstart,LaStart,L1Start,L2Start,tg,ts,1);
[FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1,L2,tg,ts,1);

fprintf('Total time of optimization %e sec \n',(t));

