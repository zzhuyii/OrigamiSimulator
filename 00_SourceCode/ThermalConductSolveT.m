%% Solve for the temperature of the system
% created: Yi Zhu 2020-04-20

function [T,indexArray]=ThermalConductSolveT(qin,ThermalMat,ThermalNode,ModelConstant,RTnode)
    A=size(RTnode);
    N=A(1);
    
    airLayer=ModelConstant{31};    
    RT=ModelConstant{23};
    
    if airLayer~=0
        T=zeros(airLayer*ThermalNode-N,1);
        Ts=RT*ones(ThermalNode+N,1);

        indexArray=[1:(airLayer*ThermalNode)];
        indexArray(RTnode)=[];

        qf=qin(indexArray);

        Kff=ThermalMat(indexArray,indexArray);
        Ksf=ThermalMat([RTnode',(airLayer*ThermalNode)+1:end],indexArray);

        T=Kff\(qf-Ksf'*Ts); 

        % the output index array has no air nodes
        indexArray=[1:ThermalNode];
        indexArray(RTnode)=[];
    else
        T=zeros(ThermalNode-N,1);
        Ts=RT*ones(N,1);
        indexArray=[1:ThermalNode];
        indexArray(RTnode)=[];
        qf=qin(indexArray);
        Kff=ThermalMat(indexArray,indexArray);
        Ksf=ThermalMat([RTnode'],indexArray);
        T=Kff\(qf-Ksf'*Ts); 
    end
end
