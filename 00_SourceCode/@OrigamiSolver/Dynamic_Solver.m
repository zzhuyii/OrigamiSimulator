function []=Dynamic_Solver(obj,Fext,supp,U0,V0)

    supp=[1,1,1,1;
          4,1,1,1;
          16,1,1,1;
          9,1,1,1;
          10,1,1,1;
          11,1,1,1;
          12,1,1,1;];
    
    dt=10^-6;
   
    U=obj.currentU;
    A=size(U);
    newNodeNum=A(1);    
    A=size(obj.sprK);
    barNum=A(1);
    
    U0=zeros(size(U));
    V0=zeros(size(U));
    
    currentAppliedForce=zeros(3*newNodeNum,1);    
    for i=1:newNodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = obj.currentAppliedForce(i,:);
    end  
   
    % for the storage of input
    Tload=zeros(newNodeNum,3);
    Sx=zeros(barNum,1);
    Ex=zeros(barNum,1);
    M=zeros(barNum,1);
    theta=zeros(barNum,1);      
    
    % Set up the mass matrix of the system
    NodalMass=obj.Dynamic_MassVector();    
    M=diag(3*newNodeNum);    
    for i=1:newNodeNum
       M(3*(i-1)+1,3*(i-1)+1)=NodalMass(i); 
       M(3*(i-1)+2,3*(i-1)+2)=NodalMass(i); 
       M(3*(i-1)+3,3*(i-1)+3)=NodalMass(i); 
    end
    M=obj.Dynamic_ModMforSupp(M,supp);
    
    
    
    for i=1:100
        
    
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
            unLoad=currentAppliedForce-Tbar-Tspr-Tcontact; 
            K=Kbar+Kspr+Kcontact; 
        else
            Tload=-(Tbar+Tspr);
            K=Kbar+Kspr;   
            unLoad=currentAppliedForce-Tbar-Tspr; 
        end

        [K,unLoad]=obj.Solver_ModKforSupp(K,supp,unLoad,nonRigidSupport,suppElastic,U);
        K=sparse(K);
        
        
    end


end