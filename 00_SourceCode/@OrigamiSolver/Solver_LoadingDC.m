%% Displacement controled method method. 
%
% this code load the structure and trace the equilibrium with a 
% displacement controled method. 'dispControler' gives the node index 
% and direction that is used as the controling displacement entries
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
%       dispControler: the controling displacment entires
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
    nodeForce,loadForce,contactForce]=Solver_LoadingDC(obj,dc)

    
    increStep=dc.increStep;
    tol=dc.tol;
    iterMax=dc.iterMax;
    lambdaBar=dc.lambdaBar;  
    
    supp=dc.supp;
    load=dc.load;
    nonRigidSupport=dc.nonRigidSupport;
    suppElastic=dc.suppElastic;    
    selectedRefDisp=dc.selectedRefDisp;
    
    newNodeNum = size(obj.newNode);
    newNodeNum = newNodeNum(1);
    currentAppliedForce=zeros(3*newNodeNum,1);    
    for i=1:newNodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = obj.currentAppliedForce(i,:);
    end   
        
    % set up the output data
    loadHis=zeros(increStep,1);
    A=size(obj.currentU);  Num=A(1);  Num2=A(2);
    UhisLoading=zeros(increStep,Num,Num2);
    U=obj.currentU;
    strainEnergyLoading=zeros(increStep,4);
    
    
    fprintf('Loading Analysis Start');
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

    selectedReferenceNodeNum=selectedRefDisp(1)*3-3+selectedRefDisp(2);
    

    for i=1:increStep
        if lastStep<=4
            lambdaBar=1*lambdaBar;
        elseif lastStep >=15
            lambdaBar=1*lambdaBar;
        end
        
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
        
        % we can deactive the GSP scalling
%         GSP=1;
%         sig=1;
        
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
                dLambda=-ur(selectedReferenceNodeNum)/up(selectedReferenceNodeNum);
            else
                dLambda=-ur(selectedReferenceNodeNum)/up(selectedReferenceNodeNum);
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
        
        loadHis(i)=norm(pload);         
        loadHis(i)=lambda;    
        
        A=size(theta);
        N=A(1);
        for j=1:N
            if obj.barType(j)==5
                strainEnergyLoading(i,2)=strainEnergyLoading(i,2)...
                    +0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                strainEnergyLoading(i,3)=strainEnergyLoading(i,3)...
                    +0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            elseif  obj.barType(j)==1
                strainEnergyLoading(i,1)=strainEnergyLoading(i,1)...
                    +0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                strainEnergyLoading(i,3)=strainEnergyLoading(i,3)...
                    +0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            else
                if obj.sprIJKL(j,1)==0
                else
                    strainEnergyLoading(i,1)=strainEnergyLoading(i,1)...
                        +0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                end
                strainEnergyLoading(i,4)=strainEnergyLoading(i,4)...
                    +0.5*obj.barArea(j)*obj.creaseE*(Ex(j))^2*obj.barLength(j);
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

