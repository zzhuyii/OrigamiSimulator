%% Modified generalized displacement controlled method. 
%
% this code load the structure and trace the equilibrium with a modified
% general displacement controled method. This method use the similar
% constrain as the cylindrical mehtod so it may experience difficulty when
% analyzing the structure self-contact. 
%
% Input:
%       panel0: the oringal panel information before meshing;
%       newNode2OldNode: mapping between the new node to the old node;
%       oldCreaseNum: total number of creases before meshing;
%       newNode: the new nodal coordinates;
%       barConnect: the two nodal coordinates of bar;
%       barType: the type of bars;
%       barArea: the area of each bar;
%       barLength: the length of each bar;
%       sprIJKL: the connectivity of each spring element;
%       sprK: the spring stiffness of each spring element;
%       sprTargetZeroStrain: the target folding stress free position;
%       creaseRef: mapping between the new bar elements to the old creases;
%       load: the applied load;
%       supportInfo: support information;
%       U: current deformation field;
%       loadConstant: loading constant to control the solver;
%       modelGeometryConstant: geometrical constant of origami;
%       modelMechanicalconstant: mechanical constat of origami;
% Output:
%       U: current diformation field;
%       UhisLoading: the history of deformation field;
%       loadHis: the applied load history;
%       strainEnergyLoading: the strain energy history;
%       nodeForce: nodal force vector;
%       loadForce: load force vector;
%       contactForce: contact force vector;
%

function [U,UhisLoading,loadHis,strainEnergyLoading,...
    nodeForce,loadForce,contactForce]...
    =Solver_LoadingMGDCM(obj,mgdcm)

    fprintf('Loading Analysis Start');

    % load the current status of the origami
    U=obj.currentU;
    newNodeNum = size(obj.newNode);
    newNodeNum = newNodeNum(1);
    currentAppliedForce=zeros(3*newNodeNum,1);    
    for i=1:newNodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = obj.currentAppliedForce(i,:);
    end   
    
    % loading controller
    increStep=mgdcm.increStep;
    tol=mgdcm.tol;
    iterMax=mgdcm.iterMax;
    lambdaBar=mgdcm.lambdaBar;  
    
    supp=mgdcm.supp;
    load=mgdcm.load;
    nonRigidSupport=mgdcm.nonRigidSupport;
    suppElastic=mgdcm.suppElastic;  
        
    % set up output data structure
    A=size(U);
    Num=A(1);
    Num2=A(2);
    loadHis=zeros(mgdcm.increStep,1);
    UhisLoading=zeros(mgdcm.increStep,Num,Num2);
    strainEnergyLoading=zeros(mgdcm.increStep,4);
    
    % Set up storage matrix for stress, strain, nodal force
    A=size(obj.sprK);
    barNum=A(1);
    mgdcm.FnodalHis=zeros(mgdcm.increStep,Num,Num2);
    mgdcm.barSxHis=zeros(mgdcm.increStep,barNum);
    mgdcm.barExHis=zeros(mgdcm.increStep,barNum);
    mgdcm.sprMHis=zeros(mgdcm.increStep,barNum);
    mgdcm.sprRotHis=zeros(mgdcm.increStep,barNum);
    
    % Assemble the load vector
    A=size(load);
    loadSize=A(1);
    loadVec=zeros(3*Num,1);
    for i=1:loadSize
        TempNodeNum=load(i,1);
        loadVec(TempNodeNum*3-2)=load(i,2);
        loadVec(TempNodeNum*3-1)=load(i,3);
        loadVec(TempNodeNum*3-0)=load(i,4);
    end
    pload=loadVec;
    up1=zeros(Num*3,increStep);
    lambda=1;
    sig=1;
    lastStep=10;

    for i=1:increStep
        if lastStep<=4
            lambdaBar=2*lambdaBar;
        elseif lastStep >=15
            lambdaBar=0.5*lambdaBar;
        end
        
        % for the storage of input
        Tload=zeros(Num,Num2);
        Sx=zeros(barNum,1);
        Ex=zeros(barNum,1);
        M=zeros(barNum,1);
        theta=zeros(barNum,1); 
        
        step=1;     
        fprintf('Icrement = %d\n',i);
        [Ex]=obj.Bar_Strain(U,obj.newNode,obj.barArea,...
            obj.barConnect,obj.barLength);
        [theta]=obj.Spr_Theta(U,obj.sprIJKL,obj.newNode);
        [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
        [M,sprKadj]=obj.Spr_Cons(obj.currentSprZeroStrain,theta,...
            obj.sprK,obj.creaseRef,obj.oldCreaseNum,...
            obj.panelInnerBarStart,obj.sprIJKL,...
            obj.newNode,U,obj.compliantCreaseOpen);
        [Kbar]=obj.Bar_GlobalStiffAssemble(U,Sx,C,obj.barArea,...
            obj.barLength,obj.barConnect,obj.newNode);
        [Kspr]=obj.Spr_GlobalStiffAssemble(U,M,obj.sprIJKL,...
            sprKadj,obj.newNode);
     
        if obj.contactOpen==1
            [Tcontact,Kcontact]=obj.Contact_AssembleForceStiffness(...
                obj.panel0,obj.newNode2OldNode,obj.newNode,...
                U,obj.ke,obj.d0edge,obj.d0center,...
                obj.centerNodeStart,obj.compliantCreaseOpen);
            K=Kbar+Kspr+Kcontact;
        else
            K=Kbar+Kspr;                
        end
      
        [K,loadVec]=obj.Solver_ModKforSupp(K,supp,loadVec,nonRigidSupport,suppElastic,U);
        up1(:,i)=K\loadVec;
        if i==1
            GSP=1;
            sig=1;
        else
            GSP=(up1(:,1)'*up1(:,1))/(up1(:,i)'*up1(:,i));
            sig=sign(up1(:,i-1)'*up1(:,i))*sig;        
        end    
        dLambda=sig*lambdaBar*sqrt(abs(GSP));
        lambda=lambda+dLambda;
        pload=pload+dLambda*loadVec;    
        dUtemp=dLambda*up1(:,i);

        for j=1:Num
            U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
        end  
        R=norm(dUtemp);
        if obj.contactOpen==1
            fprintf('    Iteration = %d, R = %e, Tcontact = %e\n',step,R,norm(Tcontact));
        else
            fprintf('    Iteration = %d, R = %e\n',step,R);
        end        
        
        step=step+1;       
        lastStep=step;

        while and(step<iterMax,R>tol)

            [Ex]=obj.Bar_Strain(U,obj.newNode,obj.barArea,...
                obj.barConnect,obj.barLength);
            [theta]=obj.Spr_Theta(U,obj.sprIJKL,obj.newNode);
            [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
            [M,sprKadj]=obj.Spr_Cons(obj.currentSprZeroStrain,theta,...
                obj.sprK,obj.creaseRef,obj.oldCreaseNum,...
                obj.panelInnerBarStart,obj.sprIJKL,obj.newNode,...
                U,obj.compliantCreaseOpen);
            [Kbar]=obj.Bar_GlobalStiffAssemble(U,Sx,C,obj.barArea,...
                obj.barLength,obj.barConnect,obj.newNode);
            [Tbar]=obj.Bar_GlobalForce(U,Sx,C,obj.barArea,...
                obj.barLength,obj.barConnect,obj.newNode);
            [Kspr]=obj.Spr_GlobalStiffAssemble(U,M,obj.sprIJKL,...
                sprKadj,obj.newNode);
            [Tspr]=obj.Spr_GlobalForce(U,M,obj.sprIJKL,...
                sprKadj,obj.newNode);            
 
            
           if obj.contactOpen==1
                [Tcontact,Kcontact]=obj.Contact_AssembleForceStiffness(...
                    obj.panel0,obj.newNode2OldNode,obj.newNode,...
                    U,obj.ke,obj.d0edge,obj.d0center,...
                    obj.centerNodeStart,obj.compliantCreaseOpen);
                Tload=-(Tbar+Tspr+Tcontact);
                unLoad=currentAppliedForce+pload-Tbar-Tspr-Tcontact; 
                K=Kbar+Kspr+Kcontact;
            else
                Tload=-(Tbar+Tspr);
                K=Kbar+Kspr;   
                unLoad=currentAppliedForce+pload-Tbar-Tspr; 
            end
            
            
            [K,unLoad]=obj.Solver_ModKforSupp(K,supp,unLoad,nonRigidSupport,suppElastic,U);

            up=K\loadVec;
            ur=K\unLoad;

            if i==1
                dLambda=-(up1(:,1)'*ur)/(up1(:,1)'*up);
            else
                dLambda=-(up1(:,i)'*ur)/(up1(:,i)'*up);
            end
            lambda=lambda+dLambda;
            pload=pload+dLambda*loadVec;

            dUtemp=dLambda*up+ur;
            for j=1:Num
                U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
            end

            R=norm(dUtemp);
            if obj.contactOpen==1
                fprintf('    Iteration = %d, R = %e, Tcontact = %e\n',step,R,norm(Tcontact));
            else
                fprintf('    Iteration = %d, R = %e\n',step,R);
            end
            
            step=step+1;  
            lastStep=step;
        end 
        % The following part is used to calculate output information
        UhisLoading(i,:,:)=U;
        
        % Store the loading history
        Tload=reshape(Tload,Num,Num2);
        mgdcm.FnodalHis(i,:,:)=Tload;
        mgdcm.barSxHis(i,:)=Sx;
        mgdcm.barExHis(i,:)=Ex;
        mgdcm.sprMHis(i,:)=M;
        mgdcm.sprRotHis(i,:)=theta;
        
        
        loadHis(i)=lambda;    
        A=size(theta);
        N=A(1);
        for j=1:N
            if obj.barType(j)==5
                strainEnergyLoading(i,2)=strainEnergyLoading(i,2)+...
                    0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                strainEnergyLoading(i,3)=strainEnergyLoading(i,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            elseif  obj.barType(j)==1
                strainEnergyLoading(i,1)=strainEnergyLoading(i,1)+...
                    0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                strainEnergyLoading(i,3)=strainEnergyLoading(i,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            else
                if obj.sprIJKL(j,1)==0
                else
                    strainEnergyLoading(i,1)=strainEnergyLoading(i,1)+...
                        0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                end
                strainEnergyLoading(i,4)=strainEnergyLoading(i,4)+...
                    0.5*obj.barArea(j)*obj.creaseE*(Ex(j))^2*obj.barLength(j);
            end          
        end    
    end
    if obj.contactOpen==1
        Tforce=Tcontact+Tbar+Tspr;
    else
        Tforce=Tbar+Tspr;
    end    
    nodeForce=zeros(size(U));
    loadForce=zeros(size(U));
    contactForce=zeros(size(U));
    for j=1:Num
        nodeForce((j),:)=(Tforce(3*j-2:3*j))';
        loadForce((j),:)=(pload(3*j-2:3*j))';
        
        if obj.contactOpen==1
            contactForce((j),:)=(Tcontact(3*j-2:3*j))';
        else
        end
        
    end
    
end

