%% Constitutive relationships for springs
% This function calculate the moment and stiffness of the rotation springs
% given the geometry, deformation, and material properties of the
% rotational springs. 
%
% Input: 
%       sprTargetZeroStrain: the target self-fold angle;
%       theta: the current folding angle of rotational springs;
%       sprK: the spring stiffness of rotational spring elements;
%       creaseRef: mapping the old crease to new bars;
%       oldCreaseNum: total number of old creases;
%       panelInnerBarStart: starting bar index of first panel bars;
%       sprIJKL: stores the connectivity of each rotational springs;
%       newNode: the nodal coordinates of each new node;
%       U: displacement field;
%       compliantCreaseOpen: index determines if compliant crease is open;
% Output:
%       M: bending moment of spring element;
%       sprKadj: spring stiffness after adjustment;
%

function  [M,sprKadj]= Spr_Cons(sprTargeZeroStrain,theta,...
    sprK,creaseRef,oldCreaseNum,panelInnerBarStart,sprIJKL,...
    newNode,U,compliantCreaseOpen)

    theta1=0*pi;
    theta2=2*pi-0*pi;

    M=zeros(size(theta));
    sprKadj=zeros(size(theta));
    A=size(theta);
    N=A(1);

    for i=panelInnerBarStart:N
        M(i)=sprK(i)*(theta(i)-sprTargeZeroStrain(i)); 
        sprKadj(i)=sprK(i);     
    end
    
    if compliantCreaseOpen==0
        
        
        for i=1:oldCreaseNum
            if sprIJKL(i,1)==0
            else
                Pi=sprIJKL(i,1);
                Pj=sprIJKL(i,2);
                Pk=sprIJKL(i,3);
                Pl=sprIJKL(i,4);
                
                nodei=newNode(Pi,:)+U(Pi,:);
                nodej=newNode(Pj,:)+U(Pj,:);
                nodek=newNode(Pk,:)+U(Pk,:);
                nodel=newNode(Pl,:)+U(Pl,:);
                
                rij=nodei-nodej;
                rkj=nodek-nodej;
                rkl=nodek-nodel;

                m=cross(rij,rkj);
                n=cross(rkj,rkl);

                if dot(m,rkl)==0
                    yita=1;
                else
                    yita=sign(dot(m,rkl));
                end
                creaseRot=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi); 
                
                if creaseRot<theta1
                    ratio=(sec(pi*(creaseRot-theta1)/2/theta1))^2;
                    sprKadj(i)=sprK(i)+sprK(i)*ratio-sprK(i);
                    M(i)=sprK(i)*(theta(i)-sprTargeZeroStrain(i))+(2*theta1*sprK(i)/pi)*tan(pi*(creaseRot-theta1)/2/theta1)-sprK(i)*creaseRot+sprK(i)*theta1;
             
                elseif creaseRot>theta2
                    ratio=(sec(pi*(creaseRot-theta2)/(4*pi-2*theta2)))^2;
                    sprKadj(i)=sprK(i)+sprK(i)*ratio-sprK(i);
                    M(i)=sprK(i)*(theta(i)-sprTargeZeroStrain(i))+(2*(2*pi-theta2)*sprK(i)/pi)*tan(pi*(creaseRot-theta2)/(4*pi-2*theta2))-sprK(i)*creaseRot+sprK(i)*theta2;
                else   
                    sprKadj(i)=sprK(i); 
                    M(i)=sprK(i)*(theta(i)-sprTargeZeroStrain(i));
                end
            end
        end
    else
        for i=1:oldCreaseNum
            if creaseRef(i,3)~=0
                topCrease=creaseRef(i,1);
                botCrease=creaseRef(i,2);
                topPanCenter=sprIJKL(topCrease,4);
                botPanCenter=sprIJKL(botCrease,4);
                creaseCenter=sprIJKL(topCrease,1);

                centerLeft=creaseRef(i,5);
                centerRight=creaseRef(i,6);

                if sprIJKL(centerLeft,2)==creaseCenter
                    pointLeft=sprIJKL(centerLeft,3);            
                else
                    pointLeft=sprIJKL(centerLeft,2);      
                end

                if sprIJKL(centerRight,2)==creaseCenter
                    pointRight=sprIJKL(centerRight,3);            
                else
                    pointRight=sprIJKL(centerRight,2);      
                end
                Pi=topPanCenter;
                Pj=pointRight;
                Pk=pointLeft;
                Pl=botPanCenter;

                nodei=newNode(Pi,:)+U(Pi,:);
                nodej=newNode(Pj,:)+U(Pj,:);
                nodek=newNode(Pk,:)+U(Pk,:);
                nodel=newNode(Pl,:)+U(Pl,:);

                rij=nodei-nodej;
                rkj=nodek-nodej;
                rkl=nodek-nodel;

                m=cross(rij,rkj);
                n=cross(rkj,rkl);

                if dot(m,rkl)==0
                    yita=1;
                else
                    yita=sign(dot(m,rkl));
                end

                creaseRot=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi); 
                HorizontalCrease=[1 2 5 6];
                DiagonalCrease=[3 4 7 8];

                
                % This part of code will only lock the horizontal rotational
                % springs of the crease

                Ktotal=0*(sprK(creaseRef(i,1))+sprK(creaseRef(i,2)) + ...
                sprK(creaseRef(i,5))+sprK(creaseRef(i,6)));
                % Ktotal=0 means local locking is disabled

                if creaseRot<theta1
                    ratio=(sec(pi*(creaseRot-theta1)/2/theta1))^2;   
                    for j=1:4                
                        index=creaseRef(i,HorizontalCrease(j));            
                        sprKadj(index)=sprK(index)+Ktotal*ratio-Ktotal;  
                        M(index)=sprK(index)*(theta(index)-sprTargeZeroStrain(index))+(2*theta1*Ktotal/pi)*tan(pi*(creaseRot-theta1)/2/theta1)-Ktotal*creaseRot+Ktotal*theta1;
                    end              
                elseif creaseRot>theta2
                    ratio=(sec(pi*(creaseRot-theta2)/(4*pi-2*theta2)))^2;
                    for j=1:4
                        index=creaseRef(i,HorizontalCrease(j));    
                        sprKadj(index)=sprK(index)+Ktotal*ratio-Ktotal;   
                        M(index)=sprK(index)*(theta(index)-sprTargeZeroStrain(index))+(2*(2*pi-theta2)*Ktotal/pi)*tan(pi*(creaseRot-theta2)/(4*pi-2*theta2))-Ktotal*creaseRot+Ktotal*theta2;
                    end 
                else
                    for j=1:4
                        index=creaseRef(i,HorizontalCrease(j));    
                        sprKadj(index)=sprK(index); 
                        M(index)=sprK(index)*(theta(index)-sprTargeZeroStrain(index));
                    end
                end

                for j=1:4
                    index=creaseRef(i,DiagonalCrease(j));    
                    sprKadj(index)=sprK(index); 
                    M(index)=sprK(index)*(theta(index)-sprTargeZeroStrain(index));            
                end
            end
        end
    end
end
