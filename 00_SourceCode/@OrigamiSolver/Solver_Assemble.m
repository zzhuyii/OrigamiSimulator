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

function [U,UhisAssemble,strainEnergyAssemble,...
            sprTargetZeroStrain,rotTargetZeroStrain]=...
            Solver_Assemble(obj,selfFold)
    
    U=obj.currentU;
    strainEnergyAssemble=zeros(selfFold.increStep,4);
    % [1] Energy stored in Creases Bending;
    % [2] Energy stored in Panel Bending;
    % [3] Energy stroed in Panel Stretching;
    % [4] Energy stored in Creases Stretching;
    A=size(U);
    Num=A(1);
    Num2=A(2);
    UhisAssemble=zeros(selfFold.increStep,Num,Num2);

    % Assemble the load vector
    A=size(obj.currentAppliedForce);
    LoadSize=A(1);
    LoadVec=zeros(3*Num,1);
    for i=1:LoadSize        
        LoadVec(3*(i-1)+1:3*i)=obj.currentAppliedForce(i,:);
    end    
    
    fprintf('Self Assemble Analysis Start');
    count=1;
    beforeSprZeroStrain=obj.currentSprZeroStrain;
    sprCurrentSprZeroStrain=beforeSprZeroStrain;
    
    sprTargetZeroStrain=obj.Mesh_CalculateZeroStrainFolding(...
        selfFold.targetRotZeroStrain);
    rotTargetZeroStrain=selfFold.targetRotZeroStrain;
    
        
    for i=1:selfFold.increStep
        % Calculate the stiffness at the begining of each Incremental
        A=size(obj.sprK);
        NumOfSpring=A(1);
        for j=1:NumOfSpring
            sprCurrentSprZeroStrain(j)=(1-i/selfFold.increStep)*beforeSprZeroStrain(j)...
                +(i/selfFold.increStep)*(sprTargetZeroStrain(j));
        end
        fprintf('Icrement = %d\n',count);
        step=1;
        R=1;
        while and(step<selfFold.iterMax,R>selfFold.tol)
            [Ex]=obj.Bar_Strain(U,obj.newNode,obj.barArea,...
                obj.barConnect,obj.barLength);
            [theta]=obj.Spr_Theta(U,obj.sprIJKL,obj.newNode);
            [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
            [M,sprKadj]=obj.Spr_Cons(sprCurrentSprZeroStrain,theta,...
                obj.sprK,obj.creaseRef,obj.oldCreaseNum,...
                obj.panelInnerBarStart,obj.sprIJKL,obj.newNode,...
                U,obj.compliantCreaseOpen);
            [Kbar]=obj.Bar_GlobalStiffAssemble(U,Sx,C,obj.barArea,...
                obj.barLength,obj.barConnect,obj.newNode);
            [Tbar]=obj.Bar_GlobalForce(U,Sx,C,obj.barArea,obj.barLength,...
                obj.barConnect,obj.newNode);
            [Kspr]=obj.Spr_GlobalStiffAssemble(U,M,obj.sprIJKL,sprKadj,...
                obj.newNode);
            [Tspr]=obj.Spr_GlobalForce(U,M,obj.sprIJKL,sprKadj,obj.newNode);

            if obj.contactOpen==1
                [Tcontact,Kcontact]=obj.Contact_AssembleForceStiffness(...
                    obj.panel0,obj.newNode2OldNode,...
                    obj.newNode,U,obj.ke,obj.d0edge,obj.d0center,...
                    obj.centerNodeStart,obj.compliantCreaseOpen);
                Tload=-(Tbar+Tspr+Tcontact);
                K=Kbar+Kspr+Kcontact;
            else
                Tload=-(Tbar+Tspr);
                K=Kbar+Kspr;                
            end 
            Tload=LoadVec+Tload;

            [K,Tload]=obj.Solver_ModKforSupp(K,selfFold.supp,Tload,...
                selfFold.nonRigidSupport,selfFold.suppElastic,U);  
            K=sparse(K);                 

            deltaU=(Tload'/K)';
            for j=1:Num
                U((j),:)=U((j),:)+deltaU(3*j-2:3*j)';
            end 

            R=norm(Tload);
            if obj.contactOpen==1
                fprintf('	Iteration = %d, R = %e, Tcontact = %e, Tspr = %e\n',step,R,norm(Tcontact),norm(Tspr));
            else
                fprintf('	Iteration = %d, R = %e, Tspr = %e\n',step,R,norm(Tspr));
            end
            step=step+1; 
        end
        A=size(theta);
        N=A(1);
        for j=1:N
            if obj.barType(j)==5
                strainEnergyAssemble(count,2)=strainEnergyAssemble(count,2)+...
                    0.5*obj.sprK(j)*(sprCurrentSprZeroStrain(j)-theta(j))^2;
                strainEnergyAssemble(count,3)=strainEnergyAssemble(count,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            elseif obj.barType(j)==1
                strainEnergyAssemble(count,1)=strainEnergyAssemble(count,1)+...
                    0.5*obj.sprK(j)*(sprCurrentSprZeroStrain(j)-theta(j))^2;
                strainEnergyAssemble(count,3)=strainEnergyAssemble(count,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            else
                if obj.sprIJKL(j,1)==0
                else
                    strainEnergyAssemble(count,1)=strainEnergyAssemble(count,1)...
                        +0.5*obj.sprK(j)*(sprCurrentSprZeroStrain(j)-theta(j))^2;
                end
                strainEnergyAssemble(count,4)=strainEnergyAssemble(count,4)+...
                    0.5*obj.barArea(j)*obj.creaseE*(Ex(j))^2*obj.barLength(j);
            end        
        end
        UhisAssemble(count,:,:)=U;
        count=count+1;
    end
end

