function [UHis]=Dynamic_Solver(obj,Fext,supp,V0,dt,Tfinal)

    % Test the solver
    supp=[1,1,1,1;
          4,1,1,1;
          16,1,1,1;
          9,1,1,1;
          10,1,1,1;
          11,1,1,1;
          12,1,1,1;];
      
    nonRigidSupport=0;
    suppElastic=0;
    
    % Set up storage    
    U=obj.currentU;
    A=size(U);
    newNodeNum=A(1);    
    A=size(obj.sprK);
    barNum=A(1);
    
    
    dt=10^-5;
    step=100000;
    TimeVec=(1:step)*dt;
    Fext=zeros(step,newNodeNum,3);
    % Apply a sin wave loading in z direction at the tips
    Fext(:,6,3)=0.000001;
    Fext(:,7,3)=0.000001;
    figure
    plot(TimeVec,Fext(:,7,3))
    
    UHis=zeros(int64(step),newNodeNum,3);
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
    NodalMass=obj.Dynamic_MassVector();    
    Mass=diag(3*newNodeNum);    
    for i=1:newNodeNum
       Mass(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       Mass(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       Mass(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    %Mass=obj.Dynamic_ModMforSupp(Mass,supp);
    %Mass=sparse(Mass);

    
    % Implement the explicit solver
    for i=1:step-1
        
        % First, assemble the stiff ness matrix and the internal forces
        [Ex]=obj.Bar_Strain(squeeze(UHis(i,:,:)),obj.newNode,...
            obj.barArea,obj.barConnect,obj.barLength);
        [theta]=obj.Spr_Theta(squeeze(UHis(i,:,:)),...
            obj.sprIJKL,obj.newNode);
        [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
        [M,sprKadj]=obj.Spr_Cons(obj.currentSprZeroStrain,theta,...
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
        alpha=0.0002;
        beta=0.0002;
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
        
        a=1;
        
    end


end