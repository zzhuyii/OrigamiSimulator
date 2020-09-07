%% Solve the folding angle with Timoshenko model
% 
% This function calcualte the rotation of a bi-material morph subjected to
% the elavating temperature using Timoshenko's model;
%
% Input: 
%       Tcrease: the nodal temperature;
%       modelTimoshenkoConstant: inputs for Timoshenko's model;
%       modelThermalConstant: inputs for heat transfer;
%       creaseWidth: crease of the width;
% Output:
%       rot: rotation angle after subjected to elavating Tcrease;
%

function rot=Thermal_Timoshenko(Tcrease,...
    modelTimoshenkoConstant,modelThermalConstant,creaseWidth)
    
    RT=modelThermalConstant{8};
    deltaT=Tcrease-RT;
    
    deltaa=modelTimoshenkoConstant{1};
    
    E1=modelTimoshenkoConstant{2};
    E2=modelTimoshenkoConstant{3};
    
    t1=modelTimoshenkoConstant{4};
    t2=modelTimoshenkoConstant{5};
    
    h=t1+t2;    
    m=t1/t2;
    n=E1/E2;
    
    dstrain=deltaT*deltaa;

    kappa=6*dstrain*(1+m)^2/h/(3*(1+m)^2+(1+m*n)*(m^2+1/m/n));
    rot=(creaseWidth*kappa);
    
end