%% Displacement controled method method. 
% this code load the structure and trace the equilibrium with a 
% displacement controled method. 

% 'DispControler' gives the node index and direction that is used as the
% controling displacement entries

function [U,Uhis,Loadhis,StrainEnergy,NodeForce,LoadForce,lockForce]...
    =NonlinearSolverLoadingDC(Panel,newNode,BarArea,BarConnect,...
    BarLength,BarType,SprIJKL,SprK,Theta0, ...
    Supp,Load,U,CreaseRef,CreaseNum,OldNode,...
    LoadConstant,ModelConstant,SuppElastic)

    IncreStep=LoadConstant(1);
    Tor=LoadConstant(2);
    iterMax=LoadConstant(3);
    LambdaBar=LoadConstant(4);
    
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
    DispControler=ModelConstant{19};

    Loadhis=zeros(IncreStep,1);

    A=size(U);
    Num=A(1);
    Num2=A(2);
    Uhis=zeros(IncreStep,Num,Num2);
    StrainEnergy=zeros(IncreStep,4);
    fprintf('Loading Analysis Start');


    % Assemble the load vector
    A=size(Load);
    LoadSize=A(1);
    LoadVec=zeros(3*Num,1);
    for i=1:LoadSize
        TempNodeNum=Load(i,1);
        LoadVec(TempNodeNum*3-2)=Load(i,2);
        LoadVec(TempNodeNum*3-1)=Load(i,3);
        LoadVec(TempNodeNum*3-0)=Load(i,4);
    end
    pload=LoadVec;
    up1=zeros(Num*3,IncreStep);
    lambda=1;
    sig=1;
    lastStep=10;

    selectedReferenceNodeNum=DispControler(1)*3-3+DispControler(2);
    

    for i=1:IncreStep
        if lastStep<=4
            LambdaBar=1*LambdaBar;
        elseif lastStep >=15
            LambdaBar=1*LambdaBar;
        end
        
        step=1;     
        fprintf('Icrement = %d\n',i);
        [Ex]=BarStrain(U,newNode,BarArea,BarConnect,BarLength);
        [Theta]=CreaseTheta(U,SprIJKL,newNode);
        [Sx,C]=BarCons(BarType,Ex,PanelE,CreaseE);
        [M,newCreaseK,Sx,C]=CreaseCons(Theta0,Theta,SprK,CreaseRef,CreaseNum,...
            PanelInerBarStart,SprIJKL,newNode,U,Sx,C,CreaseE,CompliantCreaseOpen);
        [Kbar]=BarGlobalAssemble(U,Sx,C,BarArea,BarLength,BarConnect,newNode);
        [Kspr]=CreaseGlobalAssemble(U,M,SprIJKL,newCreaseK,newNode);
     
        if LockingOpen==1
            [Tlock,Klock]=LockingAssemble(Panel,newNode,...
                U,CenterNodeStart,CreaseW,OldNode,ke,...
                d0edge,d0center,CompliantCreaseOpen);
            K=Kbar+Kspr+Klock;
        else
            K=Kbar+Kspr;                
        end
        
%         if i==1
%             figure;
%             spy(K);                  
%         end
%         if i==30
%             figure;
%             spy(K);                  
%         end
        
        [K,LoadVec]=ModKforSupp(K,Supp,LoadVec,ElasticSupportOpen,SuppElastic);
        up1(:,i)=K\LoadVec;  
        
        if i==1
            GSP=1;
            sig=1;
        else
            GSP=(up1(:,1)'*up1(:,1))/(up1(:,i)'*up1(:,i));
            sig=sign(up1(:,i-1)'*up1(:,i))*sig;        
        end
        
%         GSP=1;
%         sig=1;
        
        dLambda=sig*LambdaBar*sqrt(abs(GSP));
        lambda=lambda+dLambda;
        pload=pload+dLambda*LoadVec;    
        dUtemp=dLambda*up1(:,i);

        for j=1:Num
            U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
        end  
        R=norm(dUtemp);
        if LockingOpen==1
            fprintf('    Iteration = %d, R = %e, Tlock = %e\n',step,R,norm(Tlock));
        else
            fprintf('    Iteration = %d, R = %e\n',step,R);
        end
        step=step+1;       
        lastStep=step;

        while and(step<iterMax,R>Tor)

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
                Tload=-(Tbar+Tspr+Tlock);
                unLoad=pload-Tbar-Tspr-Tlock; 
                K=Kbar+Kspr+Klock;
            else
                Tload=-(Tbar+Tspr);
                K=Kbar+Kspr;   
                unLoad=pload-Tbar-Tspr; 
            end
            
            
            [K,unLoad]=ModKforSupp(K,Supp,unLoad,ElasticSupportOpen,SuppElastic);

            up=K\LoadVec;
            ur=K\unLoad;

            if i==1
                dLambda=-ur(selectedReferenceNodeNum)/up(selectedReferenceNodeNum);
            else
                dLambda=-ur(selectedReferenceNodeNum)/up(selectedReferenceNodeNum);
            end
            lambda=lambda+dLambda;
            pload=pload+dLambda*LoadVec;
            
            dUtemp=dLambda*up+ur;
            for j=1:Num
                U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
            end

            R=norm(dUtemp);
            if LockingOpen==1
                fprintf('    Iteration = %d, R = %e, Tlock = %e\n',step,R,norm(Tlock));
            else
                fprintf('    Iteration = %d, R = %e\n',step,R);
            end
            step=step+1;  
            lastStep=step;
        end 
        % The following part is used to calculate output information
        Uhis(i,:,:)=U;
        
        Loadhis(i)=norm(pload);         
        Loadhis(i)=lambda;    
        
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
        LoadForce((j),:)=(pload(3*j-2:3*j))';
        
        if LockingOpen==1
            lockForce((j),:)=(Tlock(3*j-2:3*j))';
        else
        end        
    end    
end

