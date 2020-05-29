%% FreqencyAnalysis

function [frequencySquared,Umode]=FrequencyAnalysis(NodalMass,U,...
            Panel,newNode,BarArea,BarConnect,BarLength,...
            BarType,SprIJKL,SprK,Supp,CreaseRef,CreaseNum,...
            OldNode,ModelConstant,SuppElastic)

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
        
    [CreaseCurrentZeroStrain]=CreaseTheta(U,SprIJKL,newNode);
    
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

    [K,Tload]=ModKforSupp(K,Supp,Tload,ElasticSupportOpen,SuppElastic,U);  
    nodeNum=size(NodalMass,1);
    M=zeros(3*nodeNum);
    for i=1:nodeNum
       M(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       M(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       M(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    M=ModMforSupp(M,Supp,ElasticSupportOpen,SuppElastic,U);
    [Umode,frequencySquared]=eig(K,M);  
end