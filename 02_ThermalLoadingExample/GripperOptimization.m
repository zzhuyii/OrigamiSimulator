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
Tmax=0.9*TgSU8;

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

[FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(Wstart,LaStart,L1Start,L2Start,tg,ts,1);
[FinalQBase,FinalQFunc,FinalTBase,FinalTFunc]=GripperFunc(W,La,L1,L2,tg,ts,1);

fprintf('Total time of optimization %e sec \n',(t));
