%% Nonlinear solver for folding (Assembling)
%
% This code folds the structure by atering the stress free angle of the
% springs incrementally and trace the equilibrium with a Newton-Raphson 
% method. This code cannot capture snapthrough during the folding. 
%
% Input:
%       panel0: the oringal panel information before meshing;
%       newNode2OldNode: mapping between the new node to the old node;
%       oldCreaseNum: total number of creases before meshing;
%       newNodeBeforeSelfFold: the nodal coordinates before self fold;
%       barConnect: the two nodal coordinates of bar;
%       barType: the type of bars;
%       barLength: the length of each bar;
%       barArea: the area of each bar;
%       sprIJKL: the connectivity of each spring element;
%       sprK: the spring stiffness of each spring element;
%       sprTargetZeroStrain: the target folding stress free position;
%       sprFoldingSequence: the folding sequence of each spring;
%       creaseRef: mapping between the new bar elements to the old creases;
%       oldCreaseNum: total number of old crease;
%       assembleConstant: constant to control the self assemble;
%       loadConstant: loading constant to control the solver;
%       modelGeometryConstant: geometrical constant of origami;
%       modelMechanicalconstant: mechanical constat of origami;
%       supportInfo: support information;
%       load: the applied load;
% Output:
%       U: current diformation field;
%       UhisLoading: the history of deformation field;
%       loadHis: the applied load history;
%       strainEnergyLoading: the strain energy history;
%       nodeForce: nodal force vector;
%       loadForce: load force vector;
%       contactForce: contact force vector;
%

function [U,UhisAssemble,strainEnergyAssemble]=Solver_Assemble(panel0,...
            newNode2OldNode,newNodeBeforeSelfFold,barConnect,barType,barLength,barArea,...
            sprIJKL,sprK,sprTargetZeroStrain,sprFoldingSequence, ...
            creaseRef,oldCreaseNum,assembleConstant,...
            modelGeometryConstant,modelMechanicalConstant,...
            supportInfo,load)

    increStep=assembleConstant(1);
    tor=assembleConstant(2);
    iterMax=assembleConstant(3);
    
    flag2D3D=modelGeometryConstant{1};
    compliantCreaseOpen=modelGeometryConstant{2};
    creaseWidthMat=modelGeometryConstant{3};
    panelInerBarStart=modelGeometryConstant{4};
    centerNodeStart=modelGeometryConstant{5};
    type1BarNum=modelGeometryConstant{6};
    
    panelE=modelMechanicalConstant{1};
    creaseE=modelMechanicalConstant{2};
    panelPoisson=modelMechanicalConstant{3};
    creasePoisson=modelMechanicalConstant{4};
    panelThicknessMat=modelMechanicalConstant{5};
    creaseThicknessMat=modelMechanicalConstant{6};
    diagonalRate=modelMechanicalConstant{7};
    panelW=modelMechanicalConstant{8};
        
    contactOpen=modelMechanicalConstant{9};
    ke=modelMechanicalConstant{10};
    d0edge=modelMechanicalConstant{11};
    d0center=modelMechanicalConstant{12};    
    
    rotationZeroStrain=modelMechanicalConstant{13};
    totalFoldingNum=modelMechanicalConstant{14};
    foldingSequence=modelMechanicalConstant{15};
    zeroStrianDistributionFactor=modelMechanicalConstant{16};
    
    supp=supportInfo{1};
    elasticSupportOpen=supportInfo{2};
    suppElastic=supportInfo{3};
    

    
    U=zeros(size(newNodeBeforeSelfFold));
    %[Theta0]=CreaseTheta(U,CreaseIJKL,newNode);
    strainEnergyAssemble=zeros(increStep*totalFoldingNum,4);
    % [1] Energy stored in Creases Bending;
    % [2] Energy stored in Panel Bending;
    % [3] Energy stroed in Panel Stretching;
    % [4] Energy stored in Creases Stretching;
    A=size(U);
    Num=A(1);
    Num2=A(2);
    UhisAssemble=zeros(increStep*totalFoldingNum,Num,Num2);

    % Assemble the load vector
    A=size(load);
    LoadSize=A(1);
    LoadVec=zeros(3*Num,1);
    for i=1:LoadSize
        TempNodeNum=(load(i,1));
        LoadVec(TempNodeNum*3-2)=load(i,2);
        LoadVec(TempNodeNum*3-1)=load(i,3);
        LoadVec(TempNodeNum*3-0)=load(i,4);
    end
    
    
    fprintf('Self Assemble Analysis Start');
    count=1;
    [ThetaNoAssemble]=Spr_Theta(U,sprIJKL,newNodeBeforeSelfFold);
    sprCurrentZeroStrain=ThetaNoAssemble;
    
    for t=1:totalFoldingNum
        A=size(sprFoldingSequence);
        NumOfSpring=A(1);
        
        for i=1:increStep
            % Calculate the stiffness at the begining of each Incremental
            for j=1:NumOfSpring
                if sprFoldingSequence(j)==t
                    sprCurrentZeroStrain(j)=(1-i/increStep)*ThetaNoAssemble(j)+(i/increStep)*(sprTargetZeroStrain(j));
                end
            end
            fprintf('Icrement = %d\n',count);
            step=1;
            R=1;
            while and(step<iterMax,R>tor)
                [Ex]=Bar_Strain(U,newNodeBeforeSelfFold,barArea,barConnect,barLength);
                [theta]=Spr_Theta(U,sprIJKL,newNodeBeforeSelfFold);
                [Sx,C]=Bar_Cons(barType,Ex,panelE,creaseE);
                [M,sprKadj]=Spr_Cons(sprCurrentZeroStrain,theta,sprK,creaseRef,oldCreaseNum,...
                    panelInerBarStart,sprIJKL,newNodeBeforeSelfFold,U,compliantCreaseOpen);
                [Kbar]=Bar_GlobalStiffAssemble(U,Sx,C,barArea,barLength,barConnect,newNodeBeforeSelfFold);
                [Tbar]=Bar_GlobalForce(U,Sx,C,barArea,barLength,barConnect,newNodeBeforeSelfFold);
                [Kspr]=Spr_GlobalStiffAssemble(U,M,sprIJKL,sprKadj,newNodeBeforeSelfFold);
                [Tspr]=Spr_GlobalForce(U,M,sprIJKL,sprKadj,newNodeBeforeSelfFold);
                
                if contactOpen==1
                    [Tcontact,Kcontact]=Contact_AssembleForceStiffness(panel0,newNode2OldNode,...
                        newNodeBeforeSelfFold,U,ke,d0edge,d0center,centerNodeStart,compliantCreaseOpen);
                    Tload=-(Tbar+Tspr+Tcontact);
                    K=Kbar+Kspr+Kcontact;
                else
                    Tload=-(Tbar+Tspr);
                    K=Kbar+Kspr;                
                end 
                Tload=LoadVec+Tload;
                
                [K,Tload]=Solver_ModKforSupp(K,supp,Tload,elasticSupportOpen,suppElastic,U);  
                K=sparse(K);                 

                deltaU=(Tload'/K)';
                for j=1:Num
                    U((j),:)=U((j),:)+deltaU(3*j-2:3*j)';
                end 

                R=norm(Tload);
                if contactOpen==1
                    fprintf('	Iteration = %d, R = %e, Tcontact = %e, Tspr = %e\n',step,R,norm(Tcontact),norm(Tspr));
                else
                    fprintf('	Iteration = %d, R = %e, Tspr = %e\n',step,R,norm(Tspr));
                end
                step=step+1; 
            end
            A=size(theta);
            N=A(1);
            for j=1:N
                if barType(j)==5
                    strainEnergyAssemble(count,2)=strainEnergyAssemble(count,2)+0.5*sprK(j)*(sprCurrentZeroStrain(j)-theta(j))^2;
                    strainEnergyAssemble(count,3)=strainEnergyAssemble(count,3)+0.5*barArea(j)*panelE*(Ex(j))^2*barLength(j);
                elseif barType(j)==1
                    strainEnergyAssemble(count,1)=strainEnergyAssemble(count,1)+0.5*sprK(j)*(sprCurrentZeroStrain(j)-theta(j))^2;
                    strainEnergyAssemble(count,3)=strainEnergyAssemble(count,3)+0.5*barArea(j)*panelE*(Ex(j))^2*barLength(j);
                else
                    if sprIJKL(j,1)==0
                    else
                        strainEnergyAssemble(count,1)=strainEnergyAssemble(count,1)+0.5*sprK(j)*(sprCurrentZeroStrain(j)-theta(j))^2;
                    end
                    strainEnergyAssemble(count,4)=strainEnergyAssemble(count,4)+0.5*barArea(j)*creaseE*(Ex(j))^2*barLength(j);
                end        
            end
            UhisAssemble(count,:,:)=U;
            count=count+1;
        end
    end
end

