%% Out put the K matrix

function K=Solver_CalcK(obj)

    K=[];
    U=obj.currentU;
    newNodeNum = size(obj.newNode);
    newNodeNum = newNodeNum(1);
    currentAppliedForce=zeros(3*newNodeNum,1);    
    for i=1:newNodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = obj.currentAppliedForce(i,:);
    end    

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
        K=Kbar+Kspr+Kcontact;
    else
        Tload=-(Tbar+Tspr);
        K=Kbar+Kspr;   
    end

    if obj.connectorOpen==1
        Tconnector=obj.Connector_GlobalForce(U, obj.newNode, obj.connectorK, obj.connectorNode);
        Kconnector=obj.Connector_Stiffness(U, obj.newNode, obj.connectorK, obj.connectorNode);      

        K=K+Kconnector;
        Tload=Tload-Tconnector;
    end
    
end
            
 