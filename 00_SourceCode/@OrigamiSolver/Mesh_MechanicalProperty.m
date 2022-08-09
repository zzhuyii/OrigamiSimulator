%% Calculate the mechanical properties of the origami
%
% Input:
%    [modelMechanicalConstant]: Stores model constant for mechanical
%       properties of origami;
%    [modelGeometryConstant]: Stores geometrical properties;
%    [oldCreaseType]: Type of creases in the old pattern;
%    [oldCreaseNum]: Total Number of creases in the old pattern;
%    [creaesRef]: relationship between old creases and new creases;
%    [barLength]: length of bars;
%    [panel0]: nodal number of each panel;
%    [barConnect]: two nodal numbers of the new bars;
%    [newNode]: the new nodal coordinates;
%
% Output:
%    [barArea]: area of each bars;
%    [sprK]: spring stiffness of each rotational springs;
%    [sprTargetZeroStrain]: target folding angle of each spring;
%


function Mesh_MechanicalProperty(obj)

% 
%     panelE=obj.panelE;
%     creaseE=obj.creaseE;
%     panelPoisson=obj.panelPoisson;
%     creasePoisson=obj.creasePoisson;
%     panelThicknessMat=obj.panelThickMat;
%     creaseThicknessMat=obj.creaseThickMat;
%     diagonalRate=obj.diagonalRate;
%     panelW=obj.panelW;
%             
%     currentRotZeroStrain=obj.currentRotZeroStrain;
%     zeroStrianDistributionFactor=obj.zeroStrianDistributionFactor;
%         
%     compliantCreaseOpen=obj.compliantCreaseOpen;
%     creaseWidthMat=obj.creaseWidthMat;
%     panelInnerBarStart=obj.panelInnerBarStart;
%     centerNodeStart=obj.centerNodeStart;
%     
    A=size(obj.barLength);
    barNum=A(1);
    A=size(obj.panel0);
    oldPanelNum=A(2);
    
    barArea=zeros(barNum,1);
    sprK=zeros(barNum,1);
    sprZeroStrain=pi*ones(barNum,1);    
    
    % when we need to consider the compliant creases
    if obj.compliantCreaseOpen==1
        % calculate the properties for bar/spr in the crease region
        for i=1:obj.oldCreaseNum
           if obj.oldCreaseType(i)>1

               creaseW=obj.creaseWidthVec(i);
               creaseThick=obj.creaseThickVec(i);
               l=obj.barLength(obj.creaseRef(i,1));
               theta0=obj.currentRotZeroStrain(i);            

               a1=CalculateA1(l,creaseThick,obj.creasePoisson);
               a2=CalculateA2(l,creaseW,creaseThick,obj.creasePoisson);
               a3=CalculateA3(l,creaseW,creaseThick,obj.creasePoisson);
               
               kSpr1=CalculateKspr1(obj.creaseE,creaseThick,creaseW);
               kSpr2=CalculateKspr2(obj.creaseE,creaseThick,creaseW,obj.diagonalRate);
           
               % spring stiffness
               sprK(obj.creaseRef(i,1))=obj.barLength(obj.creaseRef(i,1))*kSpr1;
               sprK(obj.creaseRef(i,2))=obj.barLength(obj.creaseRef(i,2))*kSpr1;               
               sprK(obj.creaseRef(i,5))=obj.barLength(obj.creaseRef(i,5))*kSpr1;
               sprK(obj.creaseRef(i,6))=obj.barLength(obj.creaseRef(i,6))*kSpr1;
               
               sprK(obj.creaseRef(i,3))=obj.barLength(obj.creaseRef(i,3))*kSpr2;
               sprK(obj.creaseRef(i,4))=obj.barLength(obj.creaseRef(i,4))*kSpr2;               
               sprK(obj.creaseRef(i,7))=obj.barLength(obj.creaseRef(i,7))*kSpr2;
               sprK(obj.creaseRef(i,8))=obj.barLength(obj.creaseRef(i,8))*kSpr2;
               
               % bar area
               barArea(obj.creaseRef(i,3)-1)=a1;
               barArea(obj.creaseRef(i,4)+1)=a1;               
               barArea(obj.creaseRef(i,7)-1)=a1;
               barArea(obj.creaseRef(i,8)+1)=a1;
               
               barArea(obj.creaseRef(i,3))=a2;
               barArea(obj.creaseRef(i,4))=a2;               
               barArea(obj.creaseRef(i,7))=a2;
               barArea(obj.creaseRef(i,8))=a2;
               
               barArea(obj.creaseRef(i,5))=a3;
               barArea(obj.creaseRef(i,6))=a3;               
            
               % asign targe zero strain of bars
               sprZeroStrain(obj.creaseRef(i,1))=...
                   (theta0-pi)*(1-obj.zeroStrianDistributionFactor)/2+pi;
               sprZeroStrain(obj.creaseRef(i,2))=...
                   (theta0-pi)*(1-obj.zeroStrianDistributionFactor)/2+pi;
               sprZeroStrain(obj.creaseRef(i,5))=...
                   (theta0-pi)*(obj.zeroStrianDistributionFactor)+pi;
               sprZeroStrain(obj.creaseRef(i,6))=...
                   (theta0-pi)*(obj.zeroStrianDistributionFactor)+pi;               
            
           end
        end
        
        % assemble the mechanical properties for panels
        panelSideNodeCount=1;
        for i=1:oldPanelNum
            A=size(obj.panel0{i});
            panelNodeNum=A(2);
            panelThick=obj.panelThickVec(i);

            lsum=0;
            area=0;
            for j=1:panelNodeNum
                barNum1=SearchCreaseNum(obj.barConnect,...
                    obj.centerNodeStart+i,panelSideNodeCount);
                node1=obj.newNode(obj.centerNodeStart+i,:);
                node2=obj.newNode(panelSideNodeCount,:);
                if j==panelNodeNum
                    barNum2=SearchCreaseNum(obj.barConnect,...
                        panelSideNodeCount,panelSideNodeCount+1-panelNodeNum);
                    node3=obj.newNode(panelSideNodeCount+1-panelNodeNum,:);
                else
                    barNum2=SearchCreaseNum(obj.barConnect,...
                        panelSideNodeCount,panelSideNodeCount+1);
                    node3=obj.newNode(panelSideNodeCount+1,:);                       
                end
                lsum=lsum+obj.barLength(barNum1);
                lsum=lsum+obj.barLength(barNum2);

                area=area+0.5*norm(cross(node1-node2,node3-node2));
                panelSideNodeCount=panelSideNodeCount+1;

            end
            aPanel=CalculateApanel(panelThick,area,obj.panelPoisson,lsum);
            kSprPanel=CalculateKsprpanel(obj.panelW,panelThick,obj.panelE);

            for j=1:panelNodeNum
                barArea(panelSideNodeCount-panelNodeNum+j-1)=...
                    barArea(panelSideNodeCount-panelNodeNum+j-1)+aPanel;
                barArea(panelSideNodeCount+obj.panelInnerBarStart-panelNodeNum+j-1)=...
                    barArea(panelSideNodeCount+obj.panelInnerBarStart-panelNodeNum+j-1)+aPanel;
                sprK(panelSideNodeCount+obj.panelInnerBarStart-panelNodeNum+j-1)=...
                    kSprPanel*obj.barLength(panelSideNodeCount+obj.panelInnerBarStart-panelNodeNum+j-1);
            end
        end 
        
    % when we donot need to consider the compliant creases
    else
        % Assemble the stiffness properties for creases
        for i=1:obj.oldCreaseNum
            if obj.oldCreaseType(i)>1
                creaseW=obj.creaseWidthVec(i);
                creaseThick=obj.creaseThickVec(i);
                l=obj.barLength(i);
                theta0=obj.currentRotZeroStrain(i);               
                kHinge=CalculateKhinge(creaseW,creaseThick,obj.creaseE);
                
                sprK(i)=l*kHinge;
                sprZeroStrain(i)=theta0;
            end             
        end 
        % Assemble the stiffness properties for panels
        for i=1:oldPanelNum
            A=size(obj.panel0{i});
            panelNodeNum=A(2);
            panelThick=obj.panelThickVec(i);
            A=obj.panel0{i};

            lsum=0;
            area=0;
            for j=1:panelNodeNum
                barNum1=SearchCreaseNum(obj.barConnect,...
                    obj.centerNodeStart+i,A(j));
                node1=obj.newNode(obj.centerNodeStart+i,:);
                node2=obj.newNode(A(j),:);
                if j==panelNodeNum
                    barNum2=SearchCreaseNum(obj.barConnect,...
                        A(1),A(j));
                    node3=obj.newNode(A(1),:);
                else
                    barNum2=SearchCreaseNum(obj.barConnect,...
                        A(j),A(j+1));
                    node3=obj.newNode(A(j+1),:);                       
                end
                lsum=lsum+obj.barLength(barNum1);
                lsum=lsum+obj.barLength(barNum2);

                area=area+0.5*norm(cross(node1-node2,node3-node2));
            end
            aPanel=CalculateApanel(panelThick,area,obj.panelPoisson,lsum);
            kSprPanel=CalculateKsprpanel(obj.panelW,panelThick,obj.panelE);

            
            for j=1:panelNodeNum
                barNum1=SearchCreaseNum(obj.barConnect,...
                    obj.centerNodeStart+i,A(j));
                if j==panelNodeNum
                    barNum2=SearchCreaseNum(obj.barConnect,...
                        A(1),A(j));
                else
                    barNum2=SearchCreaseNum(obj.barConnect,...
                        A(j),A(j+1));                      
                end
                barArea(barNum1)=barArea(barNum1)+aPanel;
                barArea(barNum2)=barArea(barNum1)+aPanel;
                sprK(barNum1)=...
                    kSprPanel*obj.barLength(barNum1);
            end
        end        
    end
    
    obj.barArea=barArea;
    obj.sprK=sprK;
    obj.currentSprZeroStrain=sprZeroStrain;
end

%% Functions used to define areas and rotational stiffness
function a1=CalculateA1(l,creaseThick,creasePoisson)
    a1=l*creaseThick/2/(1-creasePoisson^2);
end

function a2=CalculateA2(l,creaseW,creaseThick,creasePoisson)
    a2=((l^2+creaseW^2)^1.5)/4/(1+creasePoisson)*creaseThick/l/creaseW;
end

function a3=CalculateA3(l,creaseW,creaseThick,creasePoisson)
    a3=((l^2+creaseW^2)^1.5)/4/(1+creasePoisson)*creaseThick/l/creaseW;
end

function kSpr1=CalculateKspr1(creaseE,creaseThick,creaseW)
    kSpr1=(creaseE*creaseThick^3)/4/creaseW;
end

function kSpr2=CalculateKspr2(creaseE,creaseThick,creaseW,diagonalRate)
    kSpr2=diagonalRate*(creaseE*creaseThick^3)/4/creaseW;
end

function aPanel=CalculateApanel(panelThick,area,panelPoisson,lsum)
    aPanel=2*panelThick*area/(1-panelPoisson)/lsum;
end

function kSprPanel=CalculateKsprpanel(creaseW,panelThick,panelE)
    kSprPanel=panelE/12*(panelThick^3)/creaseW;
end

% This is for model with concentrated hinge
function kHinge=CalculateKhinge(creaseW,creaseThick,creaseE)
    kHinge=creaseE/12*(creaseThick^3)/creaseW;
end

% This is used to get the crease number
function creaseNum=SearchCreaseNum(newBarConnect,i,j)
    minIndex=min(i,j);
    maxIndex=max(i,j);
    A=size(newBarConnect);
    N=A(1);
    creaseNum=0;
    
    for i=1:N
        if min(newBarConnect(i,1),newBarConnect(i,2))==minIndex &&...
           max(newBarConnect(i,1),newBarConnect(i,2))==maxIndex
            creaseNum=i;
        end
    end
end