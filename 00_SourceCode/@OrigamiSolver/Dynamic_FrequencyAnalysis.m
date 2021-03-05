%% FreqencyAnalysis

function [frequencySquared,Umode]=Dynamic_FrequencyAnalysis(...
    obj,frequencyController)
        
    [Ex]=obj.Bar_Strain(obj.currentU,obj.newNode,obj.barArea,...
        obj.barConnect,obj.barLength);
    
    [Sx,C]=obj.Bar_Cons(obj.barType,Ex,obj.panelE,obj.creaseE);
    
    [Theta]=obj.Spr_Theta(obj.currentU,obj.sprIJKL,obj.newNode);
        
    [M,newCreaseK]=obj.Spr_Cons(obj.currentSprZeroStrain,Theta,...
        obj.sprK,obj.creaseRef,obj.oldCreaseNum,...
        obj.panelInnerBarStart,obj.sprIJKL,obj.newNode,...
        obj.currentU,obj.compliantCreaseOpen);
    
    [Kbar]=obj.Bar_GlobalStiffAssemble(obj.currentU,Sx,C,...
        obj.barArea,obj.barLength,obj.barConnect,obj.newNode);
    
    [Tbar]=obj.Bar_GlobalForce(obj.currentU,Sx,C,obj.barArea,...
        obj.barLength,obj.barConnect,obj.newNode);
    
    [Kspr]=obj.Spr_GlobalStiffAssemble(obj.currentU,M,...
        obj.sprIJKL,newCreaseK,obj.newNode);
    
    [Tspr]=obj.Spr_GlobalForce(obj.currentU,M,obj.sprIJKL,...
        newCreaseK,obj.newNode);

    if obj.contactOpen==1
        [Tlock,Klock]=obj.Contact_AssembleForceStiffness(...
            obj.panel0,obj.newNode2OldNode, obj.newNode,...
            obj.currentU,obj.ke,obj.d0edge,obj.d0center,...
            obj.centerNodeStart, obj.compliantCreaseOpen);
        Tload=-(Tbar+Tspr+Tlock);
        K=Kbar+Kspr+Klock;
    else
        Tload=-(Tbar+Tspr);
        K=Kbar+Kspr;                
    end 

    v=(obj.currentAppliedForce)';
    v=v(:);
    Tload=Tload+v;
    
    Supp=frequencyController.supp;
    ElasticSupportOpen=frequencyController.nonRigidSupport;
    SuppElastic=frequencyController.suppElastic;
    
    NodalMass=obj.Dynamic_MassVector();
    
    [K,Tload]=obj.Solver_ModKforSupp(K,Supp,Tload,...
        ElasticSupportOpen,SuppElastic,obj.currentU);  
    
    nodeNum=size(NodalMass,1);
    
    M=zeros(3*nodeNum);
    
    for i=1:nodeNum
       M(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       M(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       M(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    M=obj.Dynamic_ModMforSupp(M,Supp);
    [Umode,frequencySquared]=eig(K,M);  
end