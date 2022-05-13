%% Solver for folding with changing ambient temperature
%
% This solver solves the folding motion of the origami structure under the
% changing ambient temperature. It directly changes the ambient temperature
% to simulate the folding motion.
%

function [U,UhisThermal,energyHisThermal,temperatureHistory,...
            rotTargetZeroStrain,sprTargetZeroStrain]=...
            Solver_LoadingChangingTemperature(obj,thermal)    
        
    % load the crrent information of the system
    U=obj.currentU;
    q=obj.currentQ;
    Tinitial=obj.RT;   
    
    nodeNum=size(U);
    nodeNum=nodeNum(1);  
    
    [thermalMat]=obj.Thermal_AssembleConductMat(thermal,U);
%    qLoad=obj.Thermal_ConvertCreaseHeat2NodeHeat(thermal.targetCreaseHeating);
        
    % set up storage matrix
    temperatureHistory=zeros(nodeNum,thermal.thermalStep);
    UhisThermal=zeros(thermal.thermalStep,nodeNum,3);
    energyHisThermal=zeros(thermal.thermalStep,4);

    % set up storage matrix for stress/nodal force/strain
    thermal.FnodalHis=zeros(thermal.thermalStep,3);
    
    for i=1:thermal.thermalStep
        

        A=size(thermal.roomTempNode);
        Nrt=A(1);

        % Assign temperature profile
        obj.RT=(thermal.targetAmbientTemperature-Tinitial)*i/(thermal.thermalStep)+Tinitial;
        
        % Solve for the temperature profile
        [T,indexArray]=obj.Thermal_SolveTemperature(...
            q,thermalMat,nodeNum,thermal);
        
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
        %[thermalMat]=obj.Thermal_AssembleConductMat(thermal,U);

    end    
    obj.RT=thermal.targetAmbientTemperature;
end

