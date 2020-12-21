%% Assemble global stiffness and inner forces for contact simulation
%
% This code assembles the global stiffness and inner forces for contact.
% The contact stiffness and inner forces will be combined with
% contributions from springs and bars. This code is not called if the 
% contact is turned off.
%
% Input:
%       panel0: the original panel information before meshing;
%       newNode2OldNode: mapping from the new nodal number to the original 
%               nodal number before meshing;
%       newNode: the new nodal coordinates;
%       U: the current displacement field;
%       ke: contact potential scaling factor;
%       d0edge,d0center: d0 at the center and the edge of panel;
%       centerNodeStart: starting nodal number at the center of panel;
%       compliantCreaseOpen: if the compliant creases are used;
% Output:
%       Tcontact: inner force vector from panel contact;
%       Kcontact: stiffness matrix from panel contact;
%

function[Tcontact,Kcontact]=Contact_AssembleForceStiffness(obj,...
    panel0,newNode2OldNode,newNode,U,ke,d0edge,d0center,...
    centerNodeStart,compliantCreaseOpen)

    A=size(panel0);
    oldPanelNum=A(2);    
    Usize=size(U);
    pNum=Usize(1);
    Tcontact=zeros(3*pNum,1);
    Kcontact=zeros(3*pNum);
    NodeCount=1;
    % This index will be used to track the point
    
    for i=1:oldPanelNum
        B=size(panel0{i});
        Ctemp=panel0{i};
        N=B(2);  
        for j=1:N
            if compliantCreaseOpen==0  
                NodeCount=Ctemp(j);
            end
            Point=(newNode(NodeCount,:)+U(NodeCount,:))';            
            % pick one node as the point within the system that is on the
            % edge of the panels            
            NodeCount2=1;    
            % This index is used to track the three nodes of triangle
            for k=1:oldPanelNum
                B=size(panel0{k}); 
                C=panel0{k};
                Num2=B(2);
                if k==i
                    for m=1:Num2
                        NodeCount2=NodeCount2+1;
                    end 
                    % if k==i then they are the same panel there should not
                    % be any interaction within one single panel
                else
                    % if k!=i, then we start calculate the panel
                    % interaction, we need to first check which nodes are
                    % involved in this process.                    
                    for m=1:Num2                        
                        centerNode=centerNodeStart+k;
                        if compliantCreaseOpen==1
                            if m==Num2
                                T1index=NodeCount2;
                                T2index=NodeCount2-Num2+1;
                                T3index=centerNode;
                            else
                                T1index=NodeCount2;
                                T2index=NodeCount2+1;
                                T3index=centerNode;
                            end
                        else % If no compliant crease is used
                            if m==1
                                T1index=C(1);
                                T2index=C(Num2);
                                T3index=centerNode;
                            else
                                T1index=C(m);
                                T2index=C(m-1);
                                T3index=centerNode;
                            end
                        end
                        NodeCount2=NodeCount2+1;
                        T1=(newNode(T1index,:)+U(T1index,:))';
                        T2=(newNode(T2index,:)+U(T2index,:))';
                        T3=(newNode(T3index,:)+U(T3index,:))'; 
                        oldnode1=newNode2OldNode(T1index);
                        oldnode2=newNode2OldNode(T2index);
                        oldnode3=newNode2OldNode(T3index);
                        oldnode0=newNode2OldNode(NodeCount);
                        
                        check=0;
                        if oldnode0==-1
                            check=1;
                        else
                            check=1;
                            if oldnode0==oldnode1
                                check=0;
                            elseif oldnode0==oldnode2
                                check=0;
                            elseif oldnode0==oldnode3
                                check=0;
                            end
                        end 
                        
                        % IF the point is on the same crease as the
                        % triangle, panel interaction is ignored
                        if check==0                           
                        else
                            [d,zone,N1,N2]=obj.Contact_P2TDistance(Point,T1,T2,T3);                            
                            if d<=d0edge                            
                                if zone==0  
                                    [Dd2x2,Ddx]=obj.Contact_DerivativeZone0(Point,T1,T2,T3);
                                    index1=3*(NodeCount)-2;
                                    index2=3*(T1index)-2;
                                    index3=3*(T2index)-2;
                                    index4=3*(T3index)-2;
                                    tempF=zeros(12,1);
                                    tempK=zeros(12);                                    
                                    
                                    phi=pi/2*(1-d/d0edge);
                                    phi0=-pi/2/d0edge;
                                    for i1=1:12
                                        tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                        for j1=1:12
                                            tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0...
                                                +(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                                -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                        end
                                    end

                                    Tcontact(index1:index1+2)=Tcontact(index1:index1+2)+tempF(1:3);
                                    Tcontact(index2:index2+2)=Tcontact(index2:index2+2)+tempF(4:6);
                                    Tcontact(index3:index3+2)=Tcontact(index3:index3+2)+tempF(7:9);
                                    Tcontact(index4:index4+2)=Tcontact(index4:index4+2)+tempF(10:12);

                                    Kcontact(index1:index1+2,index1:index1+2)=Kcontact(index1:index1+2,index1:index1+2)+tempK(1:3,1:3);
                                    Kcontact(index1:index1+2,index2:index2+2)=Kcontact(index1:index1+2,index2:index2+2)+tempK(1:3,4:6);
                                    Kcontact(index1:index1+2,index3:index3+2)=Kcontact(index1:index1+2,index3:index3+2)+tempK(1:3,7:9);
                                    Kcontact(index1:index1+2,index4:index4+2)=Kcontact(index1:index1+2,index4:index4+2)+tempK(1:3,10:12);

                                    Kcontact(index2:index2+2,index1:index1+2)=Kcontact(index2:index2+2,index1:index1+2)+tempK(4:6,1:3);
                                    Kcontact(index2:index2+2,index2:index2+2)=Kcontact(index2:index2+2,index2:index2+2)+tempK(4:6,4:6);
                                    Kcontact(index2:index2+2,index3:index3+2)=Kcontact(index2:index2+2,index3:index3+2)+tempK(4:6,7:9);
                                    Kcontact(index2:index2+2,index4:index4+2)=Kcontact(index2:index2+2,index4:index4+2)+tempK(4:6,10:12);

                                    Kcontact(index3:index3+2,index1:index1+2)=Kcontact(index3:index3+2,index1:index1+2)+tempK(7:9,1:3);
                                    Kcontact(index3:index3+2,index2:index2+2)=Kcontact(index3:index3+2,index2:index2+2)+tempK(7:9,4:6);
                                    Kcontact(index3:index3+2,index3:index3+2)=Kcontact(index3:index3+2,index3:index3+2)+tempK(7:9,7:9);
                                    Kcontact(index3:index3+2,index4:index4+2)=Kcontact(index3:index3+2,index4:index4+2)+tempK(7:9,10:12);

                                    Kcontact(index4:index4+2,index1:index1+2)=Kcontact(index4:index4+2,index1:index1+2)+tempK(10:12,1:3);
                                    Kcontact(index4:index4+2,index2:index2+2)=Kcontact(index4:index4+2,index2:index2+2)+tempK(10:12,4:6);
                                    Kcontact(index4:index4+2,index3:index3+2)=Kcontact(index4:index4+2,index3:index3+2)+tempK(10:12,7:9);
                                    Kcontact(index4:index4+2,index4:index4+2)=Kcontact(index4:index4+2,index4:index4+2)+tempK(10:12,10:12);

                                elseif zone==1
                                    if N1==1
                                        tempNode=T1;
                                        index1=3*(T1index)-2;
                                    elseif N1==2
                                        tempNode=T2;
                                        index1=3*(T2index)-2;
                                    elseif N1==3
                                        tempNode=T3;
                                        index1=3*(T3index)-2;
                                    end
                                    index0=3*(NodeCount)-2;

                                    [Dd2x2,Ddx]=obj.Contact_DerivativeZone1(Point,tempNode,d);
                                    tempF=zeros(6,1);
                                    tempK=zeros(6);
                                    
                                    phi=pi/2*(1-d/d0edge);
                                    phi0=-pi/2/d0edge;
                                    for i1=1:6
                                        tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                        for j1=1:6
                                            tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                                -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                        end
                                    end

                                    Tcontact(index0:index0+2)=Tcontact(index0:index0+2)+tempF(1:3);
                                    Tcontact(index1:index1+2)=Tcontact(index1:index1+2)+tempF(4:6);

                                    Kcontact(index0:index0+2,index0:index0+2)=Kcontact(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                                    Kcontact(index0:index0+2,index1:index1+2)=Kcontact(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                                    Kcontact(index1:index1+2,index0:index0+2)=Kcontact(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                                    Kcontact(index1:index1+2,index1:index1+2)=Kcontact(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);                              

                                elseif zone==2

                                    index0=3*(NodeCount)-2;
                                    if N1==1
                                        tempNode1=T1;
                                        index1=3*(T1index)-2;
                                    elseif N1==2
                                        tempNode1=T2;
                                        index1=3*(T2index)-2;
                                    elseif N1==3
                                        tempNode1=T3;
                                        index1=3*(T3index)-2;
                                    end

                                    if N2==1
                                        tempNode2=T1;
                                        index2=3*(T1index)-2;
                                    elseif N2==2
                                        tempNode2=T2;
                                        index2=3*(T2index)-2;
                                    elseif N2==3
                                        tempNode2=T3;
                                        index2=3*(T3index)-2;
                                    end

                                    [Dd2x2,Ddx]=obj.Contact_DerivativeZone2(Point,tempNode1,tempNode2);

                                    tempF=zeros(9,1);
                                    tempK=zeros(9);

                                    phi=pi/2*(1-d/d0edge);
                                    phi0=-pi/2/d0edge;
                                    for i1=1:9
                                        tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                        for j1=1:9
                                            tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                                -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                        end
                                    end  

                                    Tcontact(index0:index0+2)=Tcontact(index0:index0+2)+tempF(1:3);
                                    Tcontact(index1:index1+2)=Tcontact(index1:index1+2)+tempF(4:6);
                                    Tcontact(index2:index2+2)=Tcontact(index2:index2+2)+tempF(7:9);

                                    Kcontact(index0:index0+2,index0:index0+2)=Kcontact(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                                    Kcontact(index0:index0+2,index1:index1+2)=Kcontact(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                                    Kcontact(index0:index0+2,index2:index2+2)=Kcontact(index0:index0+2,index2:index2+2)+tempK(1:3,7:9);

                                    Kcontact(index1:index1+2,index0:index0+2)=Kcontact(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                                    Kcontact(index1:index1+2,index1:index1+2)=Kcontact(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);
                                    Kcontact(index1:index1+2,index2:index2+2)=Kcontact(index1:index1+2,index2:index2+2)+tempK(4:6,7:9);

                                    Kcontact(index2:index2+2,index0:index0+2)=Kcontact(index2:index2+2,index0:index0+2)+tempK(7:9,1:3);
                                    Kcontact(index2:index2+2,index1:index1+2)=Kcontact(index2:index2+2,index1:index1+2)+tempK(7:9,4:6);
                                    Kcontact(index2:index2+2,index2:index2+2)=Kcontact(index2:index2+2,index2:index2+2)+tempK(7:9,7:9);                                
                                end
                            end
                        end
                    end  
                end
            end 
            NodeCount=NodeCount+1;
        end  
        
        % the following code calulate the cases when the center node of
        % panel is the point
        NodeCount2=1;
        index=centerNodeStart+i;
        Point=(newNode(index,:)+U(index,:))';
        for k=1:oldPanelNum
            B=size(panel0{k});   
            C=panel0{k};
            Num2=B(2);                
            if k==i
                for m=1:Num2
                    NodeCount2=NodeCount2+1;
                end  
            else
                for m=1:Num2
                    centerNode=centerNodeStart+k;                    

                    if compliantCreaseOpen==1
                        if m==Num2
                            T1index=NodeCount2;
                            T2index=NodeCount2-Num2+1;
                            T3index=centerNode;
                        else
                            T1index=NodeCount2;
                            T2index=NodeCount2+1;
                            T3index=centerNode;
                        end
                    else % If no compliant crease is used
                        if m==1
                            T1index=C(1);
                            T2index=C(Num2);
                            T3index=centerNode;
                        else
                            T1index=C(m);
                            T2index=C(m-1);
                            T3index=centerNode;
                        end
                    end
                    NodeCount2=NodeCount2+1;

                    T1=(newNode(T1index,:)+U(T1index,:))';
                    T2=(newNode(T2index,:)+U(T2index,:))';
                    T3=(newNode(T3index,:)+U(T3index,:))';                         

                    [d,zone,N1,N2]=obj.Contact_P2TDistance(Point,T1,T2,T3);  
                    if d<=d0center
                        if zone==0  
                            [Dd2x2,Ddx]=obj.Contact_DerivativeZone0(Point,T1,T2,T3);
                            index1=3*(index)-2;
                            index2=3*(T1index)-2;
                            index3=3*(T2index)-2;
                            index4=3*(T3index)-2;
                            tempF=zeros(12,1);
                            tempK=zeros(12);
                            
                            phi=pi/2*(1-d/d0center);
                            phi0=-pi/2/d0center;
                            for i1=1:12
                                tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                for j1=1:12
                                    tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                        -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                end
                            end

                            Tcontact(index1:index1+2)=Tcontact(index1:index1+2)+tempF(1:3);
                            Tcontact(index2:index2+2)=Tcontact(index2:index2+2)+tempF(4:6);
                            Tcontact(index3:index3+2)=Tcontact(index3:index3+2)+tempF(7:9);
                            Tcontact(index4:index4+2)=Tcontact(index4:index4+2)+tempF(10:12);

                            Kcontact(index1:index1+2,index1:index1+2)=Kcontact(index1:index1+2,index1:index1+2)+tempK(1:3,1:3);
                            Kcontact(index1:index1+2,index2:index2+2)=Kcontact(index1:index1+2,index2:index2+2)+tempK(1:3,4:6);
                            Kcontact(index1:index1+2,index3:index3+2)=Kcontact(index1:index1+2,index3:index3+2)+tempK(1:3,7:9);
                            Kcontact(index1:index1+2,index4:index4+2)=Kcontact(index1:index1+2,index4:index4+2)+tempK(1:3,10:12);

                            Kcontact(index2:index2+2,index1:index1+2)=Kcontact(index2:index2+2,index1:index1+2)+tempK(4:6,1:3);
                            Kcontact(index2:index2+2,index2:index2+2)=Kcontact(index2:index2+2,index2:index2+2)+tempK(4:6,4:6);
                            Kcontact(index2:index2+2,index3:index3+2)=Kcontact(index2:index2+2,index3:index3+2)+tempK(4:6,7:9);
                            Kcontact(index2:index2+2,index4:index4+2)=Kcontact(index2:index2+2,index4:index4+2)+tempK(4:6,10:12);

                            Kcontact(index3:index3+2,index1:index1+2)=Kcontact(index3:index3+2,index1:index1+2)+tempK(7:9,1:3);
                            Kcontact(index3:index3+2,index2:index2+2)=Kcontact(index3:index3+2,index2:index2+2)+tempK(7:9,4:6);
                            Kcontact(index3:index3+2,index3:index3+2)=Kcontact(index3:index3+2,index3:index3+2)+tempK(7:9,7:9);
                            Kcontact(index3:index3+2,index4:index4+2)=Kcontact(index3:index3+2,index4:index4+2)+tempK(7:9,10:12);

                            Kcontact(index4:index4+2,index1:index1+2)=Kcontact(index4:index4+2,index1:index1+2)+tempK(10:12,1:3);
                            Kcontact(index4:index4+2,index2:index2+2)=Kcontact(index4:index4+2,index2:index2+2)+tempK(10:12,4:6);
                            Kcontact(index4:index4+2,index3:index3+2)=Kcontact(index4:index4+2,index3:index3+2)+tempK(10:12,7:9);
                            Kcontact(index4:index4+2,index4:index4+2)=Kcontact(index4:index4+2,index4:index4+2)+tempK(10:12,10:12);

                        elseif zone==1
                            if N1==1
                                tempNode=T1;
                                index1=3*(T1index)-2;
                            elseif N1==2
                                tempNode=T2;
                                index1=3*(T2index)-2;
                            elseif N1==3
                                tempNode=T3;
                                index1=3*(T3index)-2;
                            end
                            index0=3*(index)-2;

                            [Dd2x2,Ddx]=obj.Contact_DerivativeZone1(Point,tempNode,d);
                            tempF=zeros(6,1);
                            tempK=zeros(6);

                            
                            phi=pi/2*(1-d/d0center);
                            phi0=-pi/2/d0center;
                            for i1=1:6
                                tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                for j1=1:6
                                    tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+...
                                        (phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                        -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                end
                            end

                            Tcontact(index0:index0+2)=Tcontact(index0:index0+2)+tempF(1:3);
                            Tcontact(index1:index1+2)=Tcontact(index1:index1+2)+tempF(4:6);

                            Kcontact(index0:index0+2,index0:index0+2)=Kcontact(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                            Kcontact(index0:index0+2,index1:index1+2)=Kcontact(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                            Kcontact(index1:index1+2,index0:index0+2)=Kcontact(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                            Kcontact(index1:index1+2,index1:index1+2)=Kcontact(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);                              

                        elseif zone==2
                            index0=3*(index)-2;
                            index1=0;
                            index2=0;
                            tempNode1=zeros(3,1);
                            tempNode2=zeros(3,1);
                            
                            if N1==1
                                tempNode1=T1;
                                index1=3*(T1index)-2;
                            elseif N1==2
                                tempNode1=T2;
                                index1=3*(T2index)-2;
                            elseif N1==3
                                tempNode1=T3;
                                index1=3*(T3index)-2;
                            end

                            if N2==1
                                tempNode2=T1;
                                index2=3*(T1index)-2;
                            elseif N2==2
                                tempNode2=T2;
                                index2=3*(T2index)-2;
                            elseif N2==3
                                tempNode2=T3;
                                index2=3*(T3index)-2;
                            end

                            [Dd2x2,Ddx]=obj.Contact_DerivativeZone2(Point,tempNode1,tempNode2);

                            tempF=zeros(9,1);
                            tempK=zeros(9);

                            phi=pi/2*(1-d/d0center);
                            phi0=-pi/2/d0center;
                            for i1=1:9
                                tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                for j1=1:9
                                    tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+...
                                        (phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                        -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                end
                            end
                            
                            Tcontact(index0:index0+2)=Tcontact(index0:index0+2)+tempF(1:3);
                            Tcontact(index1:index1+2)=Tcontact(index1:index1+2)+tempF(4:6);
                            Tcontact(index2:index2+2)=Tcontact(index2:index2+2)+tempF(7:9);

                            Kcontact(index0:index0+2,index0:index0+2)=Kcontact(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                            Kcontact(index0:index0+2,index1:index1+2)=Kcontact(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                            Kcontact(index0:index0+2,index2:index2+2)=Kcontact(index0:index0+2,index2:index2+2)+tempK(1:3,7:9);

                            Kcontact(index1:index1+2,index0:index0+2)=Kcontact(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                            Kcontact(index1:index1+2,index1:index1+2)=Kcontact(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);
                            Kcontact(index1:index1+2,index2:index2+2)=Kcontact(index1:index1+2,index2:index2+2)+tempK(4:6,7:9);

                            Kcontact(index2:index2+2,index0:index0+2)=Kcontact(index2:index2+2,index0:index0+2)+tempK(7:9,1:3);
                            Kcontact(index2:index2+2,index1:index1+2)=Kcontact(index2:index2+2,index1:index1+2)+tempK(7:9,4:6);
                            Kcontact(index2:index2+2,index2:index2+2)=Kcontact(index2:index2+2,index2:index2+2)+tempK(7:9,7:9);                                
                        end
                    end                        
                end  
            end
        end        
    end
end
