%% Newton-Raphson solver for loading
% This code is the Newton-Raphson sovler for loading the structure. The
% NR solve cannot capture snap-through behaviors of the system. But this
% algorithm is stable to captuer the "Locking" of structure.

function [U,Uhis,Loadhis,StrainEnergy,NodeForce,LoadForce,lockForce]=...
    DynamicSolver(Panel,newNode,BarArea,...
    BarConnect,BarLength,BarType,SprIJKL,SprK,Theta0, ...
    Supp,Load,U,CreaseRef,CreaseNum,OldNode,...
    LoadConstant,ModelConstant,SuppElastic,...
    NodalMass)

    % load constants
    IncreStep=LoadConstant(1);
    h=LoadConstant(2);
    iterMax=LoadConstant(3); %(Not used)
    epsilon=LoadConstant(4); %(Not used)
    gamma=LoadConstant(5);
    beta=LoadConstant(6);
    RayleighM=LoadConstant(7);
    RayleighK=LoadConstant(8);
    Rayleigh=LoadConstant(9);
    
    % model constants
    CreaseW=ModelConstant{1};
    PanelE=ModelConstant{2};
    CreaseE=ModelConstant{3};
    PanelThick=ModelConstant{4};
    CreaseThick=ModelConstant{5};
    PanelPoisson=ModelConstant{6};
    CreasePoisson=ModelConstant{7};
    Flag2D3D=ModelConstant{8};
    DiagonalRate=ModelConstant{9};
    LockingOpen=ModelConstant{10};
    ke=ModelConstant{11};
    d0edge=ModelConstant{12};
    d0center=ModelConstant{13};
    TotalFoldingNum=ModelConstant{14};
    PanelInerBarStart=ModelConstant{15};
    CenterNodeStart=ModelConstant{16};
    CompliantCreaseOpen=ModelConstant{17};
    ElasticSupportOpen=ModelConstant{18};
    
    % Initialize output data matrices and vectors
    Loadhis=zeros(IncreStep,1);
    
    % Assemble the mass matrix
    nodeNum=size(NodalMass,1);
    Mass=zeros(3*nodeNum);
    for i=1:nodeNum
       Mass(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       Mass(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       Mass(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    
    % Assemble the displacement matrix
    A=size(U);    Num=A(1);    Num2=A(2);
    Uhis=zeros(IncreStep,Num,Num2);
    Uhis(1,:,:)=U;
    UarrayHis=zeros(IncreStep,Num*3);
    UdotArrayHis=zeros(IncreStep,Num*3);  
    UdotdotArrayHis=zeros(IncreStep,Num*3);
    TforceHis=zeros(IncreStep,Num*3);
    StrainEnergy=zeros(IncreStep,4);
    fprintf('Loading Analysis Start \n');
    
    % Formulate the load vector and deltaF
    A=size(Load{1});
    LoadSize=A(1);
    LoadVec=zeros(IncreStep,3*Num);
    deltaLoadVec=zeros(IncreStep,3*Num);
    for i=1:IncreStep
        TempLoadi=Load{i};        
        for j=1:LoadSize
            TempNodeNum=(TempLoadi(j,1));
            LoadVec(i,TempNodeNum*3-2)=TempLoadi(j,2);
            LoadVec(i,TempNodeNum*3-1)=TempLoadi(j,3);
            LoadVec(i,TempNodeNum*3-0)=TempLoadi(j,4);
        end
    end
    for i=1:IncreStep-1
        deltaLoadVec(i,:)=LoadVec(i+1,:)-LoadVec(i,:);
    end

    % Assemble the load vector
    unbalance=zeros(IncreStep,3*Num);
    
    % Central Difference Method for dynamic loading
    % The central difference assume the fisrt two step to have displacement
    % that are zero vectors.
    for i=1:IncreStep-1               
        % calculate the internal forces and stiffness
        dUtemp=zeros(3*Num,1);
        
        % Calculate the first itertaion
        [Ex]=BarStrain(U,newNode,BarArea,BarConnect,BarLength);
        [Theta]=CreaseTheta(U,SprIJKL,newNode);
        [Sx,C]=BarCons(BarType,Ex,PanelE,CreaseE);
        [M,newCreaseK,Sx,C]=CreaseCons(Theta0,Theta,SprK,CreaseRef,CreaseNum,...
            PanelInerBarStart,SprIJKL,newNode,U,Sx,C,CreaseE,CompliantCreaseOpen);
        [Kbar]=BarGlobalAssemble(U,Sx,C,BarArea,BarLength,BarConnect,newNode);
        [Tbar]=BarGlobalForce(U,Sx,C,BarArea,BarLength,BarConnect,newNode);
        [Kspr]=CreaseGlobalAssemble(U,M,SprIJKL,newCreaseK,newNode);
        [Tspr]=CreaseGlobalForce(U,M,SprIJKL,newCreaseK,newNode);
        
        if LockingOpen==1
            [Tlock,Klock]=LockingAssemble(Panel,newNode,...
                U,CenterNodeStart,CreaseW,OldNode,ke,d0edge,d0center,CompliantCreaseOpen);
            TforceHis(i,:)=(Tbar+Tspr+Tlock)';
            K=Kbar+Kspr+Klock; 
        else
            K=Kbar+Kspr;   
            TforceHis(i,:)=(Tbar+Tspr)';
        end
        
        % Rayleigh Damping 
        Damp=RayleighM*Mass+RayleighK*K+Rayleigh*eye(size(Mass));
        
        % solve Acceleration at current step
        tempMass=Mass;
        unbalance(i,:)=(Damp*UdotArrayHis(i,:)'+TforceHis(i,:)'-LoadVec(i,:)')'; 
        [tempMass,unbalanceTemp]=ModKforSupp(tempMass,Supp,unbalance(i,:)',ElasticSupportOpen,SuppElastic,U);
        tempMass=sparse(tempMass);
        tempAcceleration=tempMass\unbalanceTemp;
        UdotdotArrayHis(i,:)=-tempAcceleration;
        
        % Form the system of equations needed to solve  
        rhs=Mass*(1/h/beta*UdotArrayHis(i,:)'+1/2/beta*UdotdotArrayHis(i,:)')+...
            Damp*(gamma/beta*UdotArrayHis(i,:)'-h*(1-gamma/2/beta)*UdotdotArrayHis(i,:)');
        RHS=rhs+deltaLoadVec(i,:)';
        LHS=1/beta/h/h*Mass+gamma/beta/h*Damp+K;      
        [LHS,RHS]=ModKforSupp(LHS,Supp,RHS,ElasticSupportOpen,SuppElastic,U);
        
        LHS=sparse(LHS);
        dUtemp=(LHS\RHS);
        
        fprintf('Dynamic Increment = %d,  dUtemp = %e,  \n' ,i,norm(dUtemp));

        % Update the next step
        for k=1:Num
            U((k),:)=squeeze(Uhis(i,k,:))'+dUtemp(3*k-2:3*k)';
        end
        
        % solve the acceleration for the next step
        UarrayHis(i+1,:)=UarrayHis(i,:)+dUtemp';
        UdotArrayHis(i+1,:)=(1-gamma/beta)*UdotArrayHis(i,:)+h*(1-gamma/2/beta)*UdotdotArrayHis(i,:)...
            +gamma/beta/h*dUtemp';   
        
        UdotdotArrayHis(i+1,:)=UdotdotArrayHis(i,:)+1/h/h/beta*dUtemp'...
            -1/h/beta*UdotArrayHis(i,:)-1/2/beta*UdotdotArrayHis(i,:);
              
        if LockingOpen==1
            fprintf('    After iteration: Tlock = %e\n',norm(Tlock));
        else
            fprintf('    After iteration: Tforce = %e\n',norm(TforceHis(i,:)));
        end
      
        % The following part is used to calculate output information 
        Uhis(i+1,:,:)=U;
        Loadhis(i)=norm(LoadVec);        
        
        A=size(Theta);
        N=A(1);
        for j=1:N
            if BarType(j)==5
                StrainEnergy(i,2)=StrainEnergy(i,2)+0.5*SprK(j)*(Theta0(j)-Theta(j))^2;
                StrainEnergy(i,3)=StrainEnergy(i,3)+0.5*BarArea(j)*PanelE*(Ex(j))^2*BarLength(j);
            elseif  BarType(j)==1
                StrainEnergy(i,1)=StrainEnergy(i,1)+0.5*SprK(j)*(Theta0(j)-Theta(j))^2;
                StrainEnergy(i,3)=StrainEnergy(i,3)+0.5*BarArea(j)*PanelE*(Ex(j))^2*BarLength(j);
            else
                if SprIJKL(j,1)==0
                else
                    StrainEnergy(i,1)=StrainEnergy(i,1)+0.5*SprK(j)*(Theta0(j)-Theta(j))^2;
                end
                StrainEnergy(i,4)=StrainEnergy(i,4)+0.5*BarArea(j)*CreaseE*(Ex(j))^2*BarLength(j);
            end         
        end
    end
    
    if LockingOpen==1
        Tforce=Tlock+Tbar+Tspr;
    else
        Tforce=Tbar+Tspr;
    end
    
    NodeForce=zeros(size(U));
    LoadForce=zeros(size(U));
    lockForce=zeros(size(U));
    for j=1:Num
        NodeForce((j),:)=(Tforce(3*j-2:3*j))';
        LoadForce((j),:)=(Tforce(3*j-2:3*j))';
        
        if LockingOpen==1
            lockForce((j),:)=(Tlock(3*j-2:3*j))';
        else
        end
        
    end    
end