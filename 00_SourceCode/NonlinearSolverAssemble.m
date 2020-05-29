%% Nonlinear solver for folding (Assembling)
% This code folds the structure by atering the stress free angle of the
% springs incrementally and trace the equilibrium with a Newton-Raphson 
% method. This code cannot capture snapthrough during the folding. 

function [U,Uhis,StrainEnergy]=NonlinearSolverAssemble(Panel,...
            newNode,BarArea,BarConnect,BarLength,...
            BarType,SprIJKL,SprK,SprTargetZeroStrain, ...
            Supp,CreaseRef,CreaseNum,NewFoldingSequence,...
            OldNode,AssembleConstant,ModelConstant,...
            SuppElastic,Load)

    IncreStep=AssembleConstant(1);
    Tor=AssembleConstant(2);
    iterMax=AssembleConstant(3);
    
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
    

    
    U=zeros(size(newNode));
    %[Theta0]=CreaseTheta(U,CreaseIJKL,newNode);
    StrainEnergy=zeros(IncreStep*TotalFoldingNum,4);
    % [1] Energy stored in Creases Bending;
    % [2] Energy stored in Panel Bending;
    % [3] Energy stroed in Panel Stretching;
    % [4] Energy stored in Creases Stretching;
    A=size(U);
    Num=A(1);
    Num2=A(2);
    Uhis=zeros(IncreStep*TotalFoldingNum,Num,Num2);

    % Assemble the load vector
    A=size(Load);
    LoadSize=A(1);
    LoadVec=zeros(3*Num,1);
    for i=1:LoadSize
        TempNodeNum=(Load(i,1));
        LoadVec(TempNodeNum*3-2)=Load(i,2);
        LoadVec(TempNodeNum*3-1)=Load(i,3);
        LoadVec(TempNodeNum*3-0)=Load(i,4);
    end
    
    
    fprintf('Self Assemble Analysis Start');
    count=1;
    [ThetaNoAssemble]=CreaseTheta(U,SprIJKL,newNode);
    CreaseCurrentZeroStrain=ThetaNoAssemble;
    
    for t=1:TotalFoldingNum
        A=size(NewFoldingSequence);
        NumOfSpring=A(1);
        
%         [Theta0]=CreaseTheta(U,SprIJKL,newNode);
%         TempTargetZeroStrain=SprTargetZeroStrain;
%         for i=1:NumOfSpring
%             if NewFoldingSequence(i)==t
%                 TempTargetZeroStrain(i)=SprTargetZeroStrain(i);
%             elseif NewFoldingSequence(i)<t
%                 TempTargetZeroStrain(i)=SprTargetZeroStrain(i);
%             else
%                 TempTargetZeroStrain(i)=ThetaNoAssemble(i);
%             end
%         end
        
        for i=1:IncreStep
            % Calculate the stiffness at the begining of each Incremental
            for j=1:NumOfSpring
                if NewFoldingSequence(j)==t
                    CreaseCurrentZeroStrain(j)=(1-i/IncreStep)*ThetaNoAssemble(j)+(i/IncreStep)*(SprTargetZeroStrain(j));
                end
            end
            fprintf('Icrement = %d\n',count);
            step=1;
            R=1;
            while and(step<iterMax,R>Tor)
                [Ex]=BarStrain(U,newNode,BarArea,BarConnect,BarLength);
                [Theta]=CreaseTheta(U,SprIJKL,newNode);
                [Sx,C]=BarCons(BarType,Ex,PanelE,CreaseE);
                [M,newCreaseK,Sx,C]=CreaseCons(CreaseCurrentZeroStrain,Theta,SprK,CreaseRef,CreaseNum,...
                    PanelInerBarStart,SprIJKL,newNode,U,Sx,C,CreaseE,CompliantCreaseOpen);
                [Kbar]=BarGlobalAssemble(U,Sx,C,BarArea,BarLength,BarConnect,newNode);
                [Tbar]=BarGlobalForce(U,Sx,C,BarArea,BarLength,BarConnect,newNode);
                [Kspr]=CreaseGlobalAssemble(U,M,SprIJKL,newCreaseK,newNode);
                [Tspr]=CreaseGlobalForce(U,M,SprIJKL,newCreaseK,newNode);
                
                if LockingOpen==1
                    [Tlock,Klock]=LockingAssemble(Panel,newNode,...
                        U,CenterNodeStart,CreaseW,...
                        OldNode,ke,d0edge,d0center,CompliantCreaseOpen);
                    Tload=-(Tbar+Tspr+Tlock);
                    K=Kbar+Kspr+Klock;
                else
                    Tload=-(Tbar+Tspr);
                    K=Kbar+Kspr;                
                end 
                Tload=LoadVec+Tload;
                
                [K,Tload]=ModKforSupp(K,Supp,Tload,ElasticSupportOpen,SuppElastic,U);  
                K=sparse(K);                 

                deltaU=(Tload'/K)';
                for j=1:Num
                    U((j),:)=U((j),:)+deltaU(3*j-2:3*j)';
                end 

                R=norm(Tload);
                if LockingOpen==1
                    fprintf('	Iteration = %d, R = %e, Tlock = %e, Tspr = %e\n',step,R,norm(Tlock),norm(Tspr));
                else
                    fprintf('	Iteration = %d, R = %e, Tspr = %e\n',step,R,norm(Tspr));
                end
                step=step+1; 
            end
            A=size(Theta);
            N=A(1);
            for j=1:N
                if BarType(j)==5
                    StrainEnergy(count,2)=StrainEnergy(count,2)+0.5*SprK(j)*(CreaseCurrentZeroStrain(j)-Theta(j))^2;
                    StrainEnergy(count,3)=StrainEnergy(count,3)+0.5*BarArea(j)*PanelE*(Ex(j))^2*BarLength(j);
                elseif BarType(j)==1
                    StrainEnergy(count,1)=StrainEnergy(count,1)+0.5*SprK(j)*(CreaseCurrentZeroStrain(j)-Theta(j))^2;
                    StrainEnergy(count,3)=StrainEnergy(count,3)+0.5*BarArea(j)*PanelE*(Ex(j))^2*BarLength(j);
                else
                    if SprIJKL(j,1)==0
                    else
                        StrainEnergy(count,1)=StrainEnergy(count,1)+0.5*SprK(j)*(CreaseCurrentZeroStrain(j)-Theta(j))^2;
                    end
                    StrainEnergy(count,4)=StrainEnergy(count,4)+0.5*BarArea(j)*CreaseE*(Ex(j))^2*BarLength(j);
                end        
            end
            Uhis(count,:,:)=U;
            count=count+1;
        end
    end
end

