%% Solve the folding angle with Timoshenko model
% created: Yi Zhu 2020-04-20

function rot=ThermalConductTimoshenko(Tnode,ModelConstant)
    
    RT=ModelConstant{23};
    deltaT=Tnode-RT;
    deltaa=ModelConstant{24};
    
    E1=ModelConstant{25};
    E2=ModelConstant{26};
    
    t1=ModelConstant{27};
    t2=ModelConstant{28};
    
    h=t1+t2;
    
    L=ModelConstant{1};
    
    m=t1/t2;
    n=E1/E2;
    
    dstrain=deltaT*deltaa;

    kappa=6*dstrain*(1+m)^2/h/(3*(1+m)^2+(1+m*n)*(m^2+1/m/n));
    rot=(L*kappa);
    
end