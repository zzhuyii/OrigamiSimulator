%% Constitutive relationships for springs
% This function calculate the moment and stiffness of the rotation springs
% given the geometry, deformation, and material properties of the
% rotational springs. Localized locking of folds can be considered here

% The content in this code is based on the open-access deformable orgami 
% simulator developed by K. Liu and G. H. Paulino 
% [1] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  

function  [M,newCreaseK,Sx,C]= CreaseCons(CreaseRotZeroStrain,Theta,SprK,CreaseRef,CreaseNum,...
    PanelInerBarStart,SprIJKL,newNode,U,Sx,C,CreaseE,CompliantCreaseOpen)

    theta1=0*pi;
    theta2=2*pi-0*pi;

    M=zeros(size(Theta));
    newCreaseK=zeros(size(Theta));
    A=size(Theta);
    N=A(1);

    for i=PanelInerBarStart:N
        M(i)=SprK(i)*(Theta(i)-CreaseRotZeroStrain(i)); 
        newCreaseK(i)=SprK(i);     
    end
    
    if CompliantCreaseOpen==0
        
        
        for i=1:CreaseNum
            if SprIJKL(i,1)==0
            else
                Pi=SprIJKL(i,1);
                Pj=SprIJKL(i,2);
                Pk=SprIJKL(i,3);
                Pl=SprIJKL(i,4);
                
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
                    newCreaseK(i)=SprK(i)+SprK(i)*ratio-SprK(i);
                    M(i)=SprK(i)*(Theta(i)-CreaseRotZeroStrain(i))+(2*theta1*SprK(i)/pi)*tan(pi*(creaseRot-theta1)/2/theta1)-SprK(i)*creaseRot+SprK(i)*theta1;
             
                elseif creaseRot>theta2
                    ratio=(sec(pi*(creaseRot-theta2)/(4*pi-2*theta2)))^2;
                    newCreaseK(i)=SprK(i)+SprK(i)*ratio-SprK(i);
                    M(i)=SprK(i)*(Theta(i)-CreaseRotZeroStrain(i))+(2*(2*pi-theta2)*SprK(i)/pi)*tan(pi*(creaseRot-theta2)/(4*pi-2*theta2))-SprK(i)*creaseRot+SprK(i)*theta2;
                else   
                    newCreaseK(i)=SprK(i); 
                    M(i)=SprK(i)*(Theta(i)-CreaseRotZeroStrain(i));
                end
            end
        end
    else
        for i=1:CreaseNum
            if CreaseRef(i,3)~=0
                topCrease=CreaseRef(i,1);
                botCrease=CreaseRef(i,2);
                topPanCenter=SprIJKL(topCrease,4);
                botPanCenter=SprIJKL(botCrease,4);
                creaseCenter=SprIJKL(topCrease,1);

                centerLeft=CreaseRef(i,5);
                centerRight=CreaseRef(i,6);

                if SprIJKL(centerLeft,2)==creaseCenter
                    pointLeft=SprIJKL(centerLeft,3);            
                else
                    pointLeft=SprIJKL(centerLeft,2);      
                end

                if SprIJKL(centerRight,2)==creaseCenter
                    pointRight=SprIJKL(centerRight,3);            
                else
                    pointRight=SprIJKL(centerRight,2);      
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

                Ktotal=0*(SprK(CreaseRef(i,1))+SprK(CreaseRef(i,2)) + ...
                SprK(CreaseRef(i,5))+SprK(CreaseRef(i,6)));
                % Ktotal=0 means local locking is disabled

                if creaseRot<theta1
                    ratio=(sec(pi*(creaseRot-theta1)/2/theta1))^2;   
                    for j=1:4                
                        index=CreaseRef(i,HorizontalCrease(j));            
                        newCreaseK(index)=SprK(index)+Ktotal*ratio-Ktotal;  
                        M(index)=SprK(index)*(Theta(index)-CreaseRotZeroStrain(index))+(2*theta1*Ktotal/pi)*tan(pi*(creaseRot-theta1)/2/theta1)-Ktotal*creaseRot+Ktotal*theta1;
                    end              
                elseif creaseRot>theta2
                    ratio=(sec(pi*(creaseRot-theta2)/(4*pi-2*theta2)))^2;
                    for j=1:4
                        index=CreaseRef(i,HorizontalCrease(j));    
                        newCreaseK(index)=SprK(index)+Ktotal*ratio-Ktotal;   
                        M(index)=SprK(index)*(Theta(index)-CreaseRotZeroStrain(index))+(2*(2*pi-theta2)*Ktotal/pi)*tan(pi*(creaseRot-theta2)/(4*pi-2*theta2))-Ktotal*creaseRot+Ktotal*theta2;
                    end 
                else
                    for j=1:4
                        index=CreaseRef(i,HorizontalCrease(j));    
                        newCreaseK(index)=SprK(index); 
                        M(index)=SprK(index)*(Theta(index)-CreaseRotZeroStrain(index));
                    end
                end

                for j=1:4
                    index=CreaseRef(i,DiagonalCrease(j));    
                    newCreaseK(index)=SprK(index); 
                    M(index)=SprK(index)*(Theta(index)-CreaseRotZeroStrain(index));            
                end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %           This section of code will lock both diagonal springs and horizontal ones 
    %           But this is not recommended
    
    %             Ktotal=0.125*(CreaseK(CreaseRef(i,1))+CreaseK(CreaseRef(i,2)) + ...
    %                         CreaseK(CreaseRef(i,5))+CreaseK(CreaseRef(i,6))+ ...
    %                          CreaseK(CreaseRef(i,3))+CreaseK(CreaseRef(i,4))+ CreaseK(CreaseRef(i,7))+CreaseK(CreaseRef(i,8)));
    % 
    %             if creaseRot<theta1
    %                 ratio=(sec(pi*(creaseRot-theta1)/2/theta1))^2;   
    %                 for j=1:8                
    %                     index=CreaseRef(i,j);            
    %                     newCreaseK(index)=CreaseK(index)+Ktotal*ratio-Ktotal;  
    %                     M(index)=CreaseK(index)*(Theta(index)-CreaseRotZeroStrain(index))+(2*theta1*Ktotal/pi)*tan(pi*(creaseRot-theta1)/2/theta1)-Ktotal*creaseRot+Ktotal*theta1;
    %                 end  
    %                 for j=1:8
    %                     index=CreaseRef(i,j);
    %                     Sx(index)=Sx(index)*(2*theta1/pi)*tan(pi*(creaseRot-theta1)/2/theta1);
    %                     C(index)=C(index)*ratio;                
    %                 end
    %             elseif creaseRot>theta2
    %                 ratio=(sec(pi*(creaseRot-theta2)/(4*pi-2*theta2)))^2;
    %                 for j=1:8
    %                     index=CreaseRef(i,j);    
    %                     newCreaseK(index)=CreaseK(index)+Ktotal*ratio-Ktotal;   
    %                     M(index)=CreaseK(index)*(Theta(index)-CreaseRotZeroStrain(index))+(2*(2*pi-theta2)*Ktotal/pi)*tan(pi*(creaseRot-theta2)/(4*pi-2*theta2))-Ktotal*creaseRot+Ktotal*theta2;
    %                 end   
    %                 for j=1:8
    %                     index=CreaseRef(i,j);
    %                     Sx(index)=Sx(index)*(2*theta1/pi)*tan(pi*(creaseRot-theta1)/2/theta1);
    %                     C(index)=C(index)*ratio;                
    %                 end
    %             else
    %                 for j=1:8
    %                     index=CreaseRef(i,j);    
    %                     newCreaseK(index)=CreaseK(index); 
    %                     M(index)=CreaseK(index)*(Theta(index)-CreaseRotZeroStrain(index));
    %                 end
    %             end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end
