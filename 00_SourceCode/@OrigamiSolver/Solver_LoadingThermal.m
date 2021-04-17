%% Nonlinear solver for folding (Assembling)
%
% This code folds the structure by atering the stress free angle of the
% springs incrementally and trace the equilibrium with a Newton-Raphson 
% method. This code cannot capture snapthrough during the folding. 
%

function [U,UhisThermal,energyHisThermal,temperatureHistory,...
            rotTargetZeroStrain,sprTargetZeroStrain]=...
            Solver_LoadingThermal(obj,thermal)    
        
    % load the crrent information of the system
    U=obj.currentU;
    q=obj.currentQ;
    T=obj.currentT;   
    
    nodeNum=size(U);
    nodeNum=nodeNum(1);  
    
    [thermalMat]=obj.Thermal_AssembleConductMat(thermal,U);
    qLoad=obj.Thermal_ConvertCreaseHeat2NodeHeat(thermal.targetCreaseHeating);
        
    % set up storage matrix
    temperatureHistory=zeros(nodeNum,thermal.thermalStep);
    UhisThermal=zeros(thermal.thermalStep,nodeNum,3);
    energyHisThermal=zeros(thermal.thermalStep,4);

    % set up storage matrix for stress/nodal force/strain
    thermal.FnodalHis=zeros(thermal.thermalStep,3);
    
    for i=1:thermal.thermalStep
            
        % linearly increment of input energy
        qtemp=i/thermal.thermalStep*qLoad+q;
        
        A=size(thermal.roomTempNode);
        Nrt=A(1);

        % Solve for the temperature profile
        [T,indexArray]=obj.Thermal_SolveTemperature(...
            qtemp,thermalMat,nodeNum,thermal);
        
        temperatureHistory(indexArray,i)=T(1:(nodeNum-Nrt));
        temperatureHistory(thermal.roomTempNode,i)=obj.RT*ones(Nrt,1);
        T=temperatureHistory(:,i);
                
        targetRot = obj.Thermal_UpdateRotZeroStrain(T,thermal);
        
        % set up the self folding process
        selfFold=ControllerSelfFolding;
        selfFold.supp=thermal.supp;
        selfFold.suppElastic=thermal.suppElastic;
        selfFold.nonRigidSupport=thermal.nonRigidSupport;
        selfFold.increStep=thermal.increStep;
        selfFold.tol=thermal.tol;
        selfFold.iterMax=thermal.iterMax;
        selfFold.plotOpen=0;
        selfFold.videoOpen=0;

        selfFold.targetRotZeroStrain = targetRot + obj.currentRotZeroStrain;        
        
        [U,UhisAssemble,strainEnergyAssemble,...
            sprTargetZeroStrain,rotTargetZeroStrain]=...
            Solver_Assemble(obj,selfFold);
        
        energyHisThermal(i,:)=squeeze(...
            strainEnergyAssemble(thermal.increStep,:)); 
        obj.currentU=U;   
        UhisThermal(i,:,:)=U;

        % update the thermal matrix for next step
        [thermalMat]=obj.Thermal_AssembleConductMat(thermal,U);

    end    
    obj.currentQ=qLoad+q;
end

