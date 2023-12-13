function [U,UHis]=Solver_Dynamics(obj,dynamics)

    % Input setup from loading controller
    % support information
    supp=dynamics.supp;
    % if elastic support is used
    nonRigidSupport=dynamics.nonRigidSupport;
    % setup of elastic support
    suppElastic=dynamics.suppElastic;
    % time step
    dt=dynamics.dt;
    % external loading forces
    Fext=dynamics.Fext;
    % total number of time steps
    A=size(Fext);
    step=A(1);
    nodeNum=A(2);
    % adjust the size of Fext
    Fext0=zeros(1,nodeNum,3);
    Fext=cat(1,Fext0,Fext);
    % vector of every time step
    TimeVec=(1:step)*dt;

    
    % Set up storage    
    U=obj.currentU;
    A=size(U);
    newNodeNum=A(1);    
    A=size(obj.sprK);
    barNum=A(1);
    
    dynamics.FnodalHis=zeros(step,newNodeNum,3);
    dynamics.barSxHis=zeros(step,barNum);
    dynamics.barExHis=zeros(step,barNum);
    dynamics.sprMHis=zeros(step,barNum);
    dynamics.sprRotHis=zeros(step,barNum);
    strainEnergy=zeros(step,4);
        
    UHis=zeros(step+1,newNodeNum,3);
    VHis=UHis;
    
    V0=zeros(size(obj.currentU));
    
    UHis(1,:,:)=obj.currentU;
    VHis(1,:,:)=V0;
    
    % The static load previouslly applied
    currentAppliedForce=zeros(3*newNodeNum,1);    
    for i=1:newNodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = obj.currentAppliedForce(i,:);
    end  
   
    % Set up the mass matrix of the system
    if obj.nodalMass==0
        obj.nodalMass=obj.Dynamic_MassVector();    
    end

    NodalMass=obj.nodalMass;    
    Mass=diag(3*newNodeNum);    
    for i=1:newNodeNum
       Mass(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       Mass(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       Mass(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    %Mass=obj.Dynamic_ModMforSupp(Mass,supp);
    %Mass=sparse(Mass);

    
    % Implement the explicit solver
    for i=1:step
        
        targetSprZeroStrain=obj.Mesh_CalculateZeroStrainFolding(...
        squeeze(dynamics.rotTargetAngle(i,:)));
        
        % First, assemble the stiff ness matrix and the internal forces
        [Ex]=obj.Bar_Strain(squeeze(UHis(i,:,:)),obj.newNode,...
            obj.barArea,obj.barConnect,obj.barLength);
        [theta]=obj.Spr_Theta(squeeze(UHis(i,:,:)),...
            obj.sprIJKL,obj.newNode);
        [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
        [M,sprKadj]=obj.Spr_Cons(targetSprZeroStrain,theta,...
            obj.sprK,obj.creaseRef,obj.oldCreaseNum,...
            obj.panelInnerBarStart,obj.sprIJKL,obj.newNode,...
            squeeze(UHis(i,:,:)),obj.compliantCreaseOpen);
        [Kbar]=obj.Bar_GlobalStiffAssemble(squeeze(UHis(i,:,:))...
            ,Sx,C,obj.barArea,obj.barLength,...
            obj.barConnect,obj.newNode);
        [Tbar]=obj.Bar_GlobalForce(squeeze(UHis(i,:,:)),Sx,C,obj.barArea,...
            obj.barLength,obj.barConnect,obj.newNode);
        [Kspr]=obj.Spr_GlobalStiffAssemble(squeeze(UHis(i,:,:)),M,obj.sprIJKL,...
            sprKadj,obj.newNode);
        [Tspr]=obj.Spr_GlobalForce(squeeze(UHis(i,:,:)),M,obj.sprIJKL,...
            sprKadj,obj.newNode);

        if obj.contactOpen==1
            [Tcontact,Kcontact]=obj.Contact_AssembleForceStiffness(...
                obj.panel0,obj.newNode2OldNode,obj.newNode,...
                squeeze(UHis(i,:,:)),obj.ke,obj.d0edge,obj.d0center,...
                obj.centerNodeStart,obj.compliantCreaseOpen);
            T=Tbar+Tspr+Tcontact;
            K=Kbar+Kspr+Kcontact; 
        else
            T=Tbar+Tspr;
            K=Kbar+Kspr;   
        end
        
        % Store the loading history
        dynamics.FnodalHis(i,:,:)=reshape(T,newNodeNum,3);
        dynamics.barSxHis(i,:)=Sx;
        dynamics.barExHis(i,:)=Ex;
        dynamics.sprMHis(i,:)=M;
        dynamics.sprRotHis(i,:)=theta;

        % Both the internal forces and the stiffness matrix is edited to
        % consider the support information
        [K,T]=obj.Solver_ModKforSupp(K,supp,T,...
            nonRigidSupport,suppElastic,UHis(i,:,:));
        K=sparse(K);

        [K,Fexti]=obj.Solver_ModKforSupp(K,supp,...
            reshape(squeeze(Fext(i,:,:))',[3*newNodeNum,1]),...
            nonRigidSupport,suppElastic,UHis(i,:,:));
        [K,Fexti1]=obj.Solver_ModKforSupp(K,supp,...
            reshape(squeeze(Fext(i+1,:,:))',[3*newNodeNum,1]),...
            nonRigidSupport,suppElastic,UHis(i,:,:));
        
        [K,Vhisi]=obj.Solver_ModKforSupp(K,supp,...
            reshape(squeeze(VHis(i,:,:))',[3*newNodeNum,1]),...
            nonRigidSupport,suppElastic,UHis(i,:,:));
        [K,Uhisi]=obj.Solver_ModKforSupp(K,supp,...
            reshape(squeeze(UHis(i,:,:))',[3*newNodeNum,1]),...
            nonRigidSupport,suppElastic,UHis(i,:,:));

            
        % Set up the damping matrix
        alpha=dynamics.alpha;
        beta=dynamics.beta;
        DampMat=alpha*Mass+beta*K;
        
        
        % Solve the acceleration
        UDotDot_i=Mass\(Fexti-DampMat*Vhisi-T);
        
        Kadjust=K+2/dt*DampMat+4/dt/dt*Mass;
        dP_adjust=(Fexti1-Fexti)+2*DampMat*Vhisi...
            +Mass*(4/dt*Vhisi+2*UDotDot_i);
        
        Uhisi1=Kadjust\dP_adjust+Uhisi;
        
        Vhisi1=2/dt*(Uhisi1-Uhisi)-Vhisi;
        
        UHis(i+1,:,:)=reshape(Uhisi1,[3,newNodeNum])';
        VHis(i+1,:,:)=reshape(Vhisi1,[3,newNodeNum])';
        
        if rem(i,1000)==0
            fprintf('finish solving %d step \n',i);
        end
        
        % Calculate strain energy in the system
        A=size(theta);
        N=A(1);
        for j=1:N
            if obj.barType(j)==5
                strainEnergy(i,2)=strainEnergy(i,2)+...
                    0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                strainEnergy(i,3)=strainEnergy(i,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            elseif  obj.barType(j)==1
                strainEnergy(i,1)=strainEnergy(i,1)+...
                    0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                strainEnergy(i,3)=strainEnergy(i,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            else
                if obj.sprIJKL(j,1)==0
                else
                    strainEnergy(i,1)=strainEnergy(i,1)+...
                        0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                end
                strainEnergy(i,4)=strainEnergy(i,4)+...
                    0.5*obj.barArea(j)*obj.creaseE*(Ex(j))^2*obj.barLength(j);
            end         
        end
        
    end
    
    dynamics.strainEnergyHis=strainEnergy;
    U=squeeze(UHis(step,:,:));
    dynamics.Uhis=UHis(1:step,:,:);

    obj.currentRotZeroStrain=dynamics.rotTargetAngle(step,:);
    obj.currentSprZeroStrain=targetSprZeroStrain;

end