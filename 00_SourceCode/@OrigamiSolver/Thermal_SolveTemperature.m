%% Solve for the temperature of the system
% 
% The function solves the temperature distribution within the active
% origmai given the input heating energy;
%
% Input: 
%       qin: input heating energy matrix;
%       thermalMat: conductivity matrix;
%       thermalNodeNum: total number of node (no environmental nodes);
%       modelThermalConstant: heat transfer input;
% Output:
%       T: temperature profile;
%       indexArray: index array with no environmental nodes;
%


function [T,indexArray]=Thermal_SolveTemperature(...
    obj,qin,thermalMat,thermalNodeNum,thermal)

    airLayer=obj.envLayer;  
    RT=obj.RT;
    roomTempNode=thermal.roomTempNode;
    A=size(roomTempNode);
    N=A(1);
    
    if airLayer~=0
        T=zeros(airLayer*thermalNodeNum-N,1);
        Ts=RT*ones(thermalNodeNum+N,1);

        indexArray=[1:(airLayer*thermalNodeNum)];
        indexArray(roomTempNode)=[];

        qf=qin(indexArray);

        Kff=thermalMat(indexArray,indexArray);
        Ksf=thermalMat([roomTempNode',(airLayer*thermalNodeNum)+1:end],indexArray);

        % Convert to sparse matrix to speed up linear equation solver
        Kff=sparse(Kff);
        
        T=Kff\(qf-Ksf'*Ts); 

        % the output index array has no air nodes
        indexArray=[1:thermalNodeNum];
        indexArray(roomTempNode)=[];
    else
        T=zeros(thermalNodeNum-N,1);
        Ts=RT*ones(N,1);
        indexArray=[1:thermalNodeNum];
        indexArray(roomTempNode)=[];
        qf=qin(indexArray);
        
        Kff=thermalMat(indexArray,indexArray);
        
        % Convert to sparse matrix to speed up linear equation solver
        Kff=sparse(Kff);
        
        Ksf=thermalMat([roomTempNode'],indexArray);
        T=Kff\(qf-Ksf'*Ts); 
    end
end
