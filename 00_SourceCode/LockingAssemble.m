%% Assemble global stiffness and global inner forces of locking
% This code assembles the global stiffness and inner forces for locking.
% The locking stiffness and inner forces will be combined with
% contributions from springs and bars.
% This code is not called if the locking is turned off.

% The content in this code is based on the contact simulation proposed 
% in the PRSA paper:
% [1] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' PRSA. (submitted)  

function[Tlock,Klock]=LockingAssemble(Panel,newNode,...
    U,CenterNodeStart,CreaseW,OldNode,ke,...
    d0edge,d0center,CompliantCreaseOpen)

    A=size(Panel);
    Npanel=A(2);    
    Usize=size(U);
    pNum=Usize(1);
    Tlock=zeros(3*pNum,1);
    Klock=zeros(3*pNum);
    NodeCount=1;
    % This index will be used to track the point
    
    for i=1:Npanel
        B=size(Panel{i});
        Ctemp=Panel{i};
        N=B(2);  
        for j=1:N
            if CompliantCreaseOpen==0  
                NodeCount=Ctemp(j);
            end
            Point=(newNode(NodeCount,:)+U(NodeCount,:))';            
            % pick one node as the point within the system that is on the
            % edge of the panels            
            NodeCount2=1;    
            % This index is used to track the three nodes of triangle
            for k=1:Npanel
                B=size(Panel{k}); 
                C=Panel{k};
                Num2=B(2);
                if k==i
                    for m=1:Num2
                        NodeCount2=NodeCount2+1;
                    end 
                    % if k==i then they are the same panel there should not
                    % be any interaction within one singel panel
                else
                    % if k!=i, then we start calculate the panel
                    % interaction, we need to first check which nodes are
                    % involved in this process.                    
                    for m=1:Num2                        
                        centerNode=CenterNodeStart+k;
                        if CompliantCreaseOpen==1
                            if m==Num2
                                T1=NodeCount2;
                                T2=NodeCount2-Num2+1;
                                T3=centerNode;
                            else
                                T1=NodeCount2;
                                T2=NodeCount2+1;
                                T3=centerNode;
                            end
                        else % If no compliant crease is used
                            if m==1
                                T1=C(1);
                                T2=C(Num2);
                                T3=centerNode;
                            else
                                T1=C(m);
                                T2=C(m-1);
                                T3=centerNode;
                            end
                        end
                        NodeCount2=NodeCount2+1;
                        pt1=(newNode(T1,:)+U(T1,:))';
                        pt2=(newNode(T2,:)+U(T2,:))';
                        pt3=(newNode(T3,:)+U(T3,:))'; 
                        oldnode1=OldNode(T1);
                        oldnode2=OldNode(T2);
                        oldnode3=OldNode(T3);
                        oldnode0=OldNode(NodeCount);
                        
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
                        
%                         InitalDistance1=norm(newNode(T1,:)-newNode(NodeCount,:));
%                         InitalDistance2=norm(newNode(T2,:)-newNode(NodeCount,:));
%                         InitalDistance3=norm(newNode(T3,:)-newNode(NodeCount,:));
%                         minArray=[InitalDistance1,InitalDistance2,InitalDistance3];
%                         minInitialDistacne=min(minArray);
                        
                        % IF the point is on the same crease as the
                        % triangle, panel interaction is ignored
                        if check==0                           
                        else
%                             Point
%                             pt1
%                             pt2
%                             pt3
                            [D,zone,N1,N2]=P2TDistance(Point,pt1,pt2,pt3);                            
                            if D<=d0edge                            
                                if zone==0  
                                    [Dd2x2,Ddx]=DerivativeZone0(Point,pt1,pt2,pt3);
                                    index1=3*(NodeCount)-2;
                                    index2=3*(T1)-2;
                                    index3=3*(T2)-2;
                                    index4=3*(T3)-2;
                                    tempF=zeros(12,1);
                                    tempK=zeros(12);                                    
                                    
                                    phi=pi/2*(1-D/d0edge);
                                    phi0=-pi/2/d0edge;
                                    for i1=1:12
                                        tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                        for j1=1:12
                                            tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0...
                                                +(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                                -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                        end
                                    end

                                    Tlock(index1:index1+2)=Tlock(index1:index1+2)+tempF(1:3);
                                    Tlock(index2:index2+2)=Tlock(index2:index2+2)+tempF(4:6);
                                    Tlock(index3:index3+2)=Tlock(index3:index3+2)+tempF(7:9);
                                    Tlock(index4:index4+2)=Tlock(index4:index4+2)+tempF(10:12);

                                    Klock(index1:index1+2,index1:index1+2)=Klock(index1:index1+2,index1:index1+2)+tempK(1:3,1:3);
                                    Klock(index1:index1+2,index2:index2+2)=Klock(index1:index1+2,index2:index2+2)+tempK(1:3,4:6);
                                    Klock(index1:index1+2,index3:index3+2)=Klock(index1:index1+2,index3:index3+2)+tempK(1:3,7:9);
                                    Klock(index1:index1+2,index4:index4+2)=Klock(index1:index1+2,index4:index4+2)+tempK(1:3,10:12);

                                    Klock(index2:index2+2,index1:index1+2)=Klock(index2:index2+2,index1:index1+2)+tempK(4:6,1:3);
                                    Klock(index2:index2+2,index2:index2+2)=Klock(index2:index2+2,index2:index2+2)+tempK(4:6,4:6);
                                    Klock(index2:index2+2,index3:index3+2)=Klock(index2:index2+2,index3:index3+2)+tempK(4:6,7:9);
                                    Klock(index2:index2+2,index4:index4+2)=Klock(index2:index2+2,index4:index4+2)+tempK(4:6,10:12);

                                    Klock(index3:index3+2,index1:index1+2)=Klock(index3:index3+2,index1:index1+2)+tempK(7:9,1:3);
                                    Klock(index3:index3+2,index2:index2+2)=Klock(index3:index3+2,index2:index2+2)+tempK(7:9,4:6);
                                    Klock(index3:index3+2,index3:index3+2)=Klock(index3:index3+2,index3:index3+2)+tempK(7:9,7:9);
                                    Klock(index3:index3+2,index4:index4+2)=Klock(index3:index3+2,index4:index4+2)+tempK(7:9,10:12);

                                    Klock(index4:index4+2,index1:index1+2)=Klock(index4:index4+2,index1:index1+2)+tempK(10:12,1:3);
                                    Klock(index4:index4+2,index2:index2+2)=Klock(index4:index4+2,index2:index2+2)+tempK(10:12,4:6);
                                    Klock(index4:index4+2,index3:index3+2)=Klock(index4:index4+2,index3:index3+2)+tempK(10:12,7:9);
                                    Klock(index4:index4+2,index4:index4+2)=Klock(index4:index4+2,index4:index4+2)+tempK(10:12,10:12);

                                elseif zone==1
                                    if N1==1
                                        tempNode=pt1;
                                        index1=3*(T1)-2;
                                    elseif N1==2
                                        tempNode=pt2;
                                        index1=3*(T2)-2;
                                    elseif N1==3
                                        tempNode=pt3;
                                        index1=3*(T3)-2;
                                    end
                                    index0=3*(NodeCount)-2;

                                    [Dd2x2,Ddx]=DerivativeZone1(Point,tempNode,D);
                                    tempF=zeros(6,1);
                                    tempK=zeros(6);
                                    
                                    phi=pi/2*(1-D/d0edge);
                                    phi0=-pi/2/d0edge;
                                    for i1=1:6
                                        tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                        for j1=1:6
                                            tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                                -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                        end
                                    end

                                    Tlock(index0:index0+2)=Tlock(index0:index0+2)+tempF(1:3);
                                    Tlock(index1:index1+2)=Tlock(index1:index1+2)+tempF(4:6);

                                    Klock(index0:index0+2,index0:index0+2)=Klock(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                                    Klock(index0:index0+2,index1:index1+2)=Klock(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                                    Klock(index1:index1+2,index0:index0+2)=Klock(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                                    Klock(index1:index1+2,index1:index1+2)=Klock(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);                              

                                elseif zone==2

                                    index0=3*(NodeCount)-2;
                                    if N1==1
                                        tempNode1=pt1;
                                        index1=3*(T1)-2;
                                    elseif N1==2
                                        tempNode1=pt2;
                                        index1=3*(T2)-2;
                                    elseif N1==3
                                        tempNode1=pt3;
                                        index1=3*(T3)-2;
                                    end

                                    if N2==1
                                        tempNode2=pt1;
                                        index2=3*(T1)-2;
                                    elseif N2==2
                                        tempNode2=pt2;
                                        index2=3*(T2)-2;
                                    elseif N2==3
                                        tempNode2=pt3;
                                        index2=3*(T3)-2;
                                    end

                                    [Dd2x2,Ddx]=DerivativeZone2(Point,tempNode1,tempNode2);

                                    tempF=zeros(9,1);
                                    tempK=zeros(9);

                                    phi=pi/2*(1-D/d0edge);
                                    phi0=-pi/2/d0edge;
                                    for i1=1:9
                                        tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                        for j1=1:9
                                            tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                                -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                        end
                                    end  

                                    Tlock(index0:index0+2)=Tlock(index0:index0+2)+tempF(1:3);
                                    Tlock(index1:index1+2)=Tlock(index1:index1+2)+tempF(4:6);
                                    Tlock(index2:index2+2)=Tlock(index2:index2+2)+tempF(7:9);

                                    Klock(index0:index0+2,index0:index0+2)=Klock(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                                    Klock(index0:index0+2,index1:index1+2)=Klock(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                                    Klock(index0:index0+2,index2:index2+2)=Klock(index0:index0+2,index2:index2+2)+tempK(1:3,7:9);

                                    Klock(index1:index1+2,index0:index0+2)=Klock(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                                    Klock(index1:index1+2,index1:index1+2)=Klock(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);
                                    Klock(index1:index1+2,index2:index2+2)=Klock(index1:index1+2,index2:index2+2)+tempK(4:6,7:9);

                                    Klock(index2:index2+2,index0:index0+2)=Klock(index2:index2+2,index0:index0+2)+tempK(7:9,1:3);
                                    Klock(index2:index2+2,index1:index1+2)=Klock(index2:index2+2,index1:index1+2)+tempK(7:9,4:6);
                                    Klock(index2:index2+2,index2:index2+2)=Klock(index2:index2+2,index2:index2+2)+tempK(7:9,7:9);                                
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
        index=CenterNodeStart+i;
        Point=(newNode(index,:)+U(index,:))';
        for k=1:Npanel
            B=size(Panel{k});   
            C=Panel{k};
            Num2=B(2);                
            if k==i
                for m=1:Num2
                    NodeCount2=NodeCount2+1;
                end  
            else
                for m=1:Num2
                    centerNode=CenterNodeStart+k;                    
%                     if m==Num2
%                         T1=NodeCount2;
%                         T2=NodeCount2-Num2+1;
%                         T3=centerNode;
%                     else
%                         T1=NodeCount2;
%                         T2=NodeCount2+1;
%                         T3=centerNode;
%                     end 
                    if CompliantCreaseOpen==1
                        if m==Num2
                            T1=NodeCount2;
                            T2=NodeCount2-Num2+1;
                            T3=centerNode;
                        else
                            T1=NodeCount2;
                            T2=NodeCount2+1;
                            T3=centerNode;
                        end
                    else % If no compliant crease is used
                        if m==1
                            T1=C(1);
                            T2=C(Num2);
                            T3=centerNode;
                        else
                            T1=C(m);
                            T2=C(m-1);
                            T3=centerNode;
                        end
                    end
                    NodeCount2=NodeCount2+1;

                    pt1=(newNode(T1,:)+U(T1,:))';
                    pt2=(newNode(T2,:)+U(T2,:))';
                    pt3=(newNode(T3,:)+U(T3,:))';                         

                    [D,zone,N1,N2]=P2TDistance(Point,pt1,pt2,pt3);  
                    if D<=d0center
                        if zone==0  
                            [Dd2x2,Ddx]=DerivativeZone0(Point,pt1,pt2,pt3);
                            index1=3*(index)-2;
                            index2=3*(T1)-2;
                            index3=3*(T2)-2;
                            index4=3*(T3)-2;
                            tempF=zeros(12,1);
                            tempK=zeros(12);
                            
                            phi=pi/2*(1-D/d0center);
                            phi0=-pi/2/d0center;
                            for i1=1:12
                                tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                for j1=1:12
                                    tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                        -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                end
                            end

                            Tlock(index1:index1+2)=Tlock(index1:index1+2)+tempF(1:3);
                            Tlock(index2:index2+2)=Tlock(index2:index2+2)+tempF(4:6);
                            Tlock(index3:index3+2)=Tlock(index3:index3+2)+tempF(7:9);
                            Tlock(index4:index4+2)=Tlock(index4:index4+2)+tempF(10:12);

                            Klock(index1:index1+2,index1:index1+2)=Klock(index1:index1+2,index1:index1+2)+tempK(1:3,1:3);
                            Klock(index1:index1+2,index2:index2+2)=Klock(index1:index1+2,index2:index2+2)+tempK(1:3,4:6);
                            Klock(index1:index1+2,index3:index3+2)=Klock(index1:index1+2,index3:index3+2)+tempK(1:3,7:9);
                            Klock(index1:index1+2,index4:index4+2)=Klock(index1:index1+2,index4:index4+2)+tempK(1:3,10:12);

                            Klock(index2:index2+2,index1:index1+2)=Klock(index2:index2+2,index1:index1+2)+tempK(4:6,1:3);
                            Klock(index2:index2+2,index2:index2+2)=Klock(index2:index2+2,index2:index2+2)+tempK(4:6,4:6);
                            Klock(index2:index2+2,index3:index3+2)=Klock(index2:index2+2,index3:index3+2)+tempK(4:6,7:9);
                            Klock(index2:index2+2,index4:index4+2)=Klock(index2:index2+2,index4:index4+2)+tempK(4:6,10:12);

                            Klock(index3:index3+2,index1:index1+2)=Klock(index3:index3+2,index1:index1+2)+tempK(7:9,1:3);
                            Klock(index3:index3+2,index2:index2+2)=Klock(index3:index3+2,index2:index2+2)+tempK(7:9,4:6);
                            Klock(index3:index3+2,index3:index3+2)=Klock(index3:index3+2,index3:index3+2)+tempK(7:9,7:9);
                            Klock(index3:index3+2,index4:index4+2)=Klock(index3:index3+2,index4:index4+2)+tempK(7:9,10:12);

                            Klock(index4:index4+2,index1:index1+2)=Klock(index4:index4+2,index1:index1+2)+tempK(10:12,1:3);
                            Klock(index4:index4+2,index2:index2+2)=Klock(index4:index4+2,index2:index2+2)+tempK(10:12,4:6);
                            Klock(index4:index4+2,index3:index3+2)=Klock(index4:index4+2,index3:index3+2)+tempK(10:12,7:9);
                            Klock(index4:index4+2,index4:index4+2)=Klock(index4:index4+2,index4:index4+2)+tempK(10:12,10:12);

                        elseif zone==1
                            if N1==1
                                tempNode=pt1;
                                index1=3*(T1)-2;
                            elseif N1==2
                                tempNode=pt2;
                                index1=3*(T2)-2;
                            elseif N1==3
                                tempNode=pt3;
                                index1=3*(T3)-2;
                            end
                            index0=3*(index)-2;

                            [Dd2x2,Ddx]=DerivativeZone1(Point,tempNode,D);
                            tempF=zeros(6,1);
                            tempK=zeros(6);

                            
                            phi=pi/2*(1-D/d0center);
                            phi0=-pi/2/d0center;
                            for i1=1:6
                                tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                for j1=1:6
                                    tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+...
                                        (phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                        -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                end
                            end

                            Tlock(index0:index0+2)=Tlock(index0:index0+2)+tempF(1:3);
                            Tlock(index1:index1+2)=Tlock(index1:index1+2)+tempF(4:6);

                            Klock(index0:index0+2,index0:index0+2)=Klock(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                            Klock(index0:index0+2,index1:index1+2)=Klock(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                            Klock(index1:index1+2,index0:index0+2)=Klock(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                            Klock(index1:index1+2,index1:index1+2)=Klock(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);                              

                        elseif zone==2
                            index0=3*(index)-2;
                            index1=0;
                            index2=0;
                            tempNode1=zeros(3,1);
                            tempNode2=zeros(3,1);
                            
                            if N1==1
                                tempNode1=pt1;
                                index1=3*(T1)-2;
                            elseif N1==2
                                tempNode1=pt2;
                                index1=3*(T2)-2;
                            elseif N1==3
                                tempNode1=pt3;
                                index1=3*(T3)-2;
                            end

                            if N2==1
                                tempNode2=pt1;
                                index2=3*(T1)-2;
                            elseif N2==2
                                tempNode2=pt2;
                                index2=3*(T2)-2;
                            elseif N2==3
                                tempNode2=pt3;
                                index2=3*(T3)-2;
                            end

                            [Dd2x2,Ddx]=DerivativeZone2(Point,tempNode1,tempNode2);

                            tempF=zeros(9,1);
                            tempK=zeros(9);

                            phi=pi/2*(1-D/d0center);
                            phi0=-pi/2/d0center;
                            for i1=1:9
                                tempF(i1)=(tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1))*ke;
                                for j1=1:9
                                    tempK(i1,j1)=((tan(phi)-phi)*Dd2x2(i1,j1)*phi0+...
                                        (phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                                        -(phi0^2)*Ddx(i1)*Ddx(j1))*ke;
                                end
                            end
                            
                            Tlock(index0:index0+2)=Tlock(index0:index0+2)+tempF(1:3);
                            Tlock(index1:index1+2)=Tlock(index1:index1+2)+tempF(4:6);
                            Tlock(index2:index2+2)=Tlock(index2:index2+2)+tempF(7:9);

                            Klock(index0:index0+2,index0:index0+2)=Klock(index0:index0+2,index0:index0+2)+tempK(1:3,1:3);
                            Klock(index0:index0+2,index1:index1+2)=Klock(index0:index0+2,index1:index1+2)+tempK(1:3,4:6);
                            Klock(index0:index0+2,index2:index2+2)=Klock(index0:index0+2,index2:index2+2)+tempK(1:3,7:9);

                            Klock(index1:index1+2,index0:index0+2)=Klock(index1:index1+2,index0:index0+2)+tempK(4:6,1:3);
                            Klock(index1:index1+2,index1:index1+2)=Klock(index1:index1+2,index1:index1+2)+tempK(4:6,4:6);
                            Klock(index1:index1+2,index2:index2+2)=Klock(index1:index1+2,index2:index2+2)+tempK(4:6,7:9);

                            Klock(index2:index2+2,index0:index0+2)=Klock(index2:index2+2,index0:index0+2)+tempK(7:9,1:3);
                            Klock(index2:index2+2,index1:index1+2)=Klock(index2:index2+2,index1:index1+2)+tempK(7:9,4:6);
                            Klock(index2:index2+2,index2:index2+2)=Klock(index2:index2+2,index2:index2+2)+tempK(7:9,7:9);                                
                        end
                    end                        
                end  
            end
        end        
    end
end
