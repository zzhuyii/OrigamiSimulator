%% Solver for dynamic electro-thermal folding 

function [U,UhisThermal,energyHisThermal,temperatureHistory,...
            rotTargetZeroStrain,sprTargetZeroStrain]=...
            Solver_DynamicsThermal(obj,dynamicThermal)    
        
    % load the crrent information of the system
    U=obj.currentU;

    supp=dynamicThermal.supp;
    nonRigidSupport=dynamicThermal.nonRigidSupport;
    suppElastic=dynamicThermal.suppElastic;
    dt=dynamicThermal.dt;
    
    nodeNum=size(U);
    nodeNum=nodeNum(1);  

    qLoadHis=obj.Thermal_ConvertCreaseHeatHis2NodeHeatHis(dynamicThermal.targetCreaseHeatingHis);
        
    % set up storage matrix
    temperatureHistory=zeros(nodeNum,dynamicThermal.step);
    energyHisThermal=zeros(dynamicThermal.step,4);

    % set up storage matrix for stress/nodal force/strain
    dynamicThermal.FnodalHis=zeros(dynamicThermal.step,3);
    
    % number of layers used to represent submerged environment
    airLayer=obj.envLayer;  

    % node number of room temperature node
    roomTempNode=dynamicThermal.roomTempNode;

    % Temperature of each step
    T=zeros((airLayer+1)*(nodeNum),1);

    % Set up storage for mechanical properties
    UhisThermal=zeros(dynamicThermal.step,nodeNum,3);
    VHis=UhisThermal;
    
    % initial velocity
    V0=zeros(size(obj.currentU));
    % deformation and velocity storage
    UhisThermal(1,:,:)=obj.currentU;
    VHis(1,:,:)=V0;
    
    % The static load previouslly applied
    currentAppliedForce=zeros(3*nodeNum,1);    
    for i=1:nodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = obj.currentAppliedForce(i,:);
    end  
   
    % Set up the mass matrix of the system
    NodalMass=obj.Dynamic_MassVector();    
    Mass=diag(3*nodeNum);    
    for i=1:nodeNum
       Mass(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       Mass(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       Mass(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    %Mass=obj.Dynamic_ModMforSupp(Mass,supp);
    %Mass=sparse(Mass);


    %% each step within the incremental solver
    for i=1:dynamicThermal.step
            
        % First we are solving the heat conduction problem
        % input heating at this dynamic loading
        qtemp=qLoadHis(:,i);

        % Assemble conductivity and capacity matrix
        [thermalMat]=obj.Thermal_AssembleConductMat(dynamicThermal,U);
        [capacityMat]=obj.Thermal_AssembleCapacityMat(dynamicThermal,U);
    
        % this represents node at the boundary
        thermalNode=size(thermalMat);
        terhmalNode=thermalNode(1);
        BCindexArray=[roomTempNode',((airLayer*nodeNum)+1):terhmalNode];

        % Convert to sparse matrix to speed up linear equation solver
        Kmodified=thermalMat;       
        Kmodified(BCindexArray,:)=0;
        Kmodified(:,BCindexArray)=0;            
        Kmodified(BCindexArray,BCindexArray)=thermalMat(BCindexArray,BCindexArray);

        matSize=size(Kmodified);
        matSize=matSize(2);

        % solve the temperature of next step            
        Cinv=1./diag(capacityMat);
        Cinv=sparse(diag(Cinv));
        % Cff is diagonal 

        A=eye(matSize)+2/3*dynamicThermal.dt*Cinv*Kmodified;
        B=eye(matSize)-1/3*dynamicThermal.dt*Cinv*Kmodified;

        C=dynamicThermal.dt*Cinv*(qtemp);
        T=A\(B*T+C);      


        % the output index array has no air nodes
        indexArray=[1:nodeNum];
        indexArray(roomTempNode)=[];

        A=size(dynamicThermal.roomTempNode);
        Nrt=A(1);

        % store the temperature solution
        temperatureHistory(indexArray,i)=T(indexArray)+obj.RT;
        temperatureHistory(dynamicThermal.roomTempNode,i)=obj.RT*ones(Nrt,1);
        

        % Now that we have a solution of temperature
        % We find the deformation of next step
        targetRot = obj.Thermal_UpdateRotZeroStrain(T,dynamicThermal) + obj.currentRotZeroStrain;
        targetSprZeroStrain=obj.Mesh_CalculateZeroStrainFolding(targetRot);
        
        % First, assemble the stiff ness matrix and the internal forces
        [Ex]=obj.Bar_Strain(squeeze(UhisThermal(i,:,:)),obj.newNode,...
            obj.barArea,obj.barConnect,obj.barLength);
        [theta]=obj.Spr_Theta(squeeze(UhisThermal(i,:,:)),...
            obj.sprIJKL,obj.newNode);
        [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
        [M,sprKadj]=obj.Spr_Cons(targetSprZeroStrain,theta,...
            obj.sprK,obj.creaseRef,obj.oldCreaseNum,...
            obj.panelInnerBarStart,obj.sprIJKL,obj.newNode,...
            squeeze(UhisThermal(i,:,:)),obj.compliantCreaseOpen);
        [Kbar]=obj.Bar_GlobalStiffAssemble(squeeze(UhisThermal(i,:,:))...
            ,Sx,C,obj.barArea,obj.barLength,...
            obj.barConnect,obj.newNode);
        [Tbar]=obj.Bar_GlobalForce(squeeze(UhisThermal(i,:,:)),Sx,C,obj.barArea,...
            obj.barLength,obj.barConnect,obj.newNode);
        [Kspr]=obj.Spr_GlobalStiffAssemble(squeeze(UhisThermal(i,:,:)),M,obj.sprIJKL,...
            sprKadj,obj.newNode);
        [Tspr]=obj.Spr_GlobalForce(squeeze(UhisThermal(i,:,:)),M,obj.sprIJKL,...
            sprKadj,obj.newNode);

        if obj.contactOpen==1
            [Tcontact,Kcontact]=obj.Contact_AssembleForceStiffness(...
                obj.panel0,obj.newNode2OldNode,obj.newNode,...
                squeeze(UhisThermal(i,:,:)),obj.ke,obj.d0edge,obj.d0center,...
                obj.centerNodeStart,obj.compliantCreaseOpen);
            Ttotal=Tbar+Tspr+Tcontact;
            K=Kbar+Kspr+Kcontact; 
        else
            Ttotal=Tbar+Tspr;
            K=Kbar+Kspr;   
        end
        
        % Store the loading history
        dynamics.FnodalHis(i,:,:)=reshape(Ttotal,nodeNum,3);
        dynamics.barSxHis(i,:)=Sx;
        dynamics.barExHis(i,:)=Ex;
        dynamics.sprMHis(i,:)=M;
        dynamics.sprRotHis(i,:)=theta;

        % Both the internal forces and the stiffness matrix is edited to
        % consider the support information
        [K,Ttotal]=obj.Solver_ModKforSupp(K,supp,Ttotal,...
            nonRigidSupport,suppElastic,UhisThermal(i,:,:));
        K=sparse(K);

        % [K,Fexti]=obj.Solver_ModKforSupp(K,supp,...
        %     reshape(squeeze(Fext(i,:,:))',[3*nodeNum,1]),...
        %     nonRigidSupport,suppElastic,UHis(i,:,:));
        % [K,Fexti1]=obj.Solver_ModKforSupp(K,supp,...
        %     reshape(squeeze(Fext(i+1,:,:))',[3*nodeNum,1]),...
        %     nonRigidSupport,suppElastic,UHis(i,:,:));
        
        [K,Vhisi]=obj.Solver_ModKforSupp(K,supp,...
            reshape(squeeze(VHis(i,:,:))',[3*nodeNum,1]),...
            nonRigidSupport,suppElastic,UhisThermal(i,:,:));
        [K,Uhisi]=obj.Solver_ModKforSupp(K,supp,...
            reshape(squeeze(UhisThermal(i,:,:))',[3*nodeNum,1]),...
            nonRigidSupport,suppElastic,UhisThermal(i,:,:));

            
        % Set up the damping matrix
        alpha=dynamicThermal.alpha;
        beta=dynamicThermal.beta;
        DampMat=alpha*Mass+beta*K;
        
        
        % Solve the acceleration
        % UDotDot_i=Mass\(Fexti-DampMat*Vhisi-T);
         % No external force for this solver

        UDotDot_i=Mass\(-DampMat*Vhisi-Ttotal);
        
        Kadjust=K+2/dt*DampMat+4/dt/dt*Mass;

        % dP_adjust=(Fexti1-Fexti)+2*DampMat*Vhisi...
        %     +Mass*(4/dt*Vhisi+2*UDotDot_i);
        % No external force for this solver

        dP_adjust=+2*DampMat*Vhisi...
            +Mass*(4/dt*Vhisi+2*UDotDot_i);
        
        Uhisi1=Kadjust\dP_adjust+Uhisi;
        
        Vhisi1=2/dt*(Uhisi1-Uhisi)-Vhisi;
        
        UhisThermal(i+1,:,:)=reshape(Uhisi1,[3,nodeNum])';
        VHis(i+1,:,:)=reshape(Vhisi1,[3,nodeNum])';

        % update the deformation
        U=squeeze(UhisThermal(i+1,:,:));
        
        if rem(i,1000)==0
            fprintf('finish solving %d step \n',i);
        end
        
        % Calculate strain energy in the system
        A=size(theta);
        N=A(1);
        for j=1:N
            if obj.barType(j)==5
                energyHisThermal(i,2)=energyHisThermal(i,2)+...
                    0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                energyHisThermal(i,3)=energyHisThermal(i,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            elseif  obj.barType(j)==1
                energyHisThermal(i,1)=energyHisThermal(i,1)+...
                    0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                energyHisThermal(i,3)=energyHisThermal(i,3)+...
                    0.5*obj.barArea(j)*obj.panelE*(Ex(j))^2*obj.barLength(j);
            else
                if obj.sprIJKL(j,1)==0
                else
                    energyHisThermal(i,1)=energyHisThermal(i,1)+...
                        0.5*obj.sprK(j)*(obj.currentSprZeroStrain(j)-theta(j))^2;
                end
                energyHisThermal(i,4)=energyHisThermal(i,4)+...
                    0.5*obj.barArea(j)*obj.creaseE*(Ex(j))^2*obj.barLength(j);
            end         
        end
        

 

    end    
    obj.currentQ=qtemp;
    dynamicThermal.strainEnergyHis=energyHisThermal;
    dynamicThermal.Uhis=UhisThermal;
    
end

