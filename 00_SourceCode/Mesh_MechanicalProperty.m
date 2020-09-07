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
%    [sprFoldingSequence]: sequence of folding of each spring;
%


function[barArea,sprK,sprTargetZeroStrain,sprFoldingSequence]...
    =Mesh_MechanicalProperty(modelMechanicalConstant,...
    modelGeometryConstant,oldCreaseType,...
    oldCreaseNum,creaseRef,barLength,panel0,...
    barConnect,newNode)


    panelE=modelMechanicalConstant{1};
    creaseE=modelMechanicalConstant{2};
    panelPoisson=modelMechanicalConstant{3};
    creasePoisson=modelMechanicalConstant{4};
    panelThicknessMat=modelMechanicalConstant{5};
    creaseThicknessMat=modelMechanicalConstant{6};
    diagonalRate=modelMechanicalConstant{7};
    panelW=modelMechanicalConstant{8};
    
    rotationZeroStrain=modelMechanicalConstant{13};
    totalFoldingNum=modelMechanicalConstant{14};
    foldingSequence=modelMechanicalConstant{15};
    zeroStrianDistributionFactor=modelMechanicalConstant{16};
    
        
    compliantCreaseOpen=modelGeometryConstant{2};
    creaseWidthMat=modelGeometryConstant{3};
    panelInerBarStart=modelGeometryConstant{4};
    centerNodeStart=modelGeometryConstant{5};
    
    A=size(barLength);
    barNum=A(1);
    A=size(panel0);
    oldPanelNum=A(2);
    
    
    barArea=zeros(barNum,1);
    sprK=zeros(barNum,1);
    sprTargetZeroStrain=pi*ones(barNum,1);
    sprFoldingSequence=zeros(barNum,1);
    
    
    % when we need to consider the compliant creases
    if compliantCreaseOpen==1
        % calculate the properties for bar/spr in the crease region
        for i=1:oldCreaseNum
           if oldCreaseType(i)>1
               creaseW=creaseWidthMat(i);
               creaseThick=creaseThicknessMat(i);
               l=barLength(creaseRef(i,1));
               theta0=rotationZeroStrain(i);               
               creaseFoldingSequence=foldingSequence(i);
                              
               a1=CalculateA1(l,creaseThick,creasePoisson);
               a2=CalculateA2(l,creaseW,creaseThick,creasePoisson);
               a3=CalculateA3(l,creaseW,creaseThick,creasePoisson);
               
               kSpr1=CalculateKspr1(creaseE,creaseThick,creaseW);
               kSpr2=CalculateKspr2(creaseE,creaseThick,creaseW,diagonalRate);
               
               % spring stiffness
               sprK(creaseRef(i,1))=barLength(creaseRef(i,1))*kSpr1;
               sprK(creaseRef(i,2))=barLength(creaseRef(i,2))*kSpr1;               
               sprK(creaseRef(i,5))=barLength(creaseRef(i,5))*kSpr1;
               sprK(creaseRef(i,6))=barLength(creaseRef(i,6))*kSpr1;
               
               sprK(creaseRef(i,3))=barLength(creaseRef(i,3))*kSpr2;
               sprK(creaseRef(i,4))=barLength(creaseRef(i,4))*kSpr2;               
               sprK(creaseRef(i,7))=barLength(creaseRef(i,7))*kSpr2;
               sprK(creaseRef(i,8))=barLength(creaseRef(i,8))*kSpr2;
               
               % bar area
               barArea(creaseRef(i,3)-1)=a1;
               barArea(creaseRef(i,4)+1)=a1;               
               barArea(creaseRef(i,7)-1)=a1;
               barArea(creaseRef(i,8)+1)=a1;
               
               barArea(creaseRef(i,3))=a2;
               barArea(creaseRef(i,4))=a2;               
               barArea(creaseRef(i,7))=a2;
               barArea(creaseRef(i,8))=a2;
               
               barArea(creaseRef(i,5))=a3;
               barArea(creaseRef(i,6))=a3;               
            
               % asign targe zero strain of bars
               sprTargetZeroStrain(creaseRef(i,1))=...
                   (theta0-pi)*(1-zeroStrianDistributionFactor)/2+pi;
               sprTargetZeroStrain(creaseRef(i,2))=...
                   (theta0-pi)*(1-zeroStrianDistributionFactor)/2+pi;
               sprTargetZeroStrain(creaseRef(i,5))=...
                   (theta0-pi)*(zeroStrianDistributionFactor)+pi;
               sprTargetZeroStrain(creaseRef(i,6))=...
                   (theta0-pi)*(zeroStrianDistributionFactor)+pi;
               
               % asign foldingSequence:
               sprFoldingSequence(creaseRef(i,1))=creaseFoldingSequence;
               sprFoldingSequence(creaseRef(i,2))=creaseFoldingSequence;
               sprFoldingSequence(creaseRef(i,5))=creaseFoldingSequence;
               sprFoldingSequence(creaseRef(i,6))=creaseFoldingSequence;               
           end
        end
        
        % assemble the mechanical properties for panels
        panelSideNodeCount=1;
        for i=1:oldPanelNum
            A=size(panel0{i});
            panelNodeNum=A(2);
            panelThick=panelThicknessMat(i);

            lsum=0;
            area=0;
            for j=1:panelNodeNum
                barNum1=SearchCreaseNum(barConnect,...
                    centerNodeStart+i,panelSideNodeCount);
                node1=newNode(centerNodeStart+i,:);
                node2=newNode(panelSideNodeCount,:);
                if j==panelNodeNum
                    barNum2=SearchCreaseNum(barConnect,...
                        panelSideNodeCount,panelSideNodeCount+1-panelNodeNum);
                    node3=newNode(panelSideNodeCount+1-panelNodeNum,:);
                else
                    barNum2=SearchCreaseNum(barConnect,...
                        panelSideNodeCount,panelSideNodeCount+1);
                    node3=newNode(panelSideNodeCount+1,:);                       
                end
                lsum=lsum+barLength(barNum1);
                lsum=lsum+barLength(barNum2);

                area=area+0.5*norm(cross(node1-node2,node3-node2));
                panelSideNodeCount=panelSideNodeCount+1;

            end
            aPanel=CalculateApanel(panelThick,area,panelPoisson,lsum);
            kSprPanel=CalculateKsprpanel(panelW,panelThick,panelE);

            for j=1:panelNodeNum
                barArea(panelSideNodeCount-panelNodeNum+j-1)=...
                    barArea(panelSideNodeCount-panelNodeNum+j-1)+aPanel;
                barArea(panelSideNodeCount+panelInerBarStart-panelNodeNum+j-1)=...
                    barArea(panelSideNodeCount+panelInerBarStart-panelNodeNum+j-1)+aPanel;
                sprK(panelSideNodeCount+panelInerBarStart-panelNodeNum+j-1)=...
                    kSprPanel*barLength(panelSideNodeCount+panelInerBarStart-panelNodeNum+j-1);
            end
        end 
        
    % when we donot need to consider the compliant creases
    else
        % Assemble the stiffness properties for creases
        for i=1:oldCreaseNum
            if oldCreaseType(i)>1
                creaseW=creaseWidthMat(i);
                creaseThick=creaseThicknessMat(i);
                l=barLength(i);
                theta0=rotationZeroStrain(i);               
                creaseFoldingSequence=foldingSequence(i);
                kHinge=CalculateKhinge(creaseW,creaseThick,creaseE);
                
                sprK(i)=l*kHinge;
                sprTargetZeroStrain(i)=theta0;
                sprFoldingSequence(i)=creaseFoldingSequence;
            end             
        end 
        % Assemble the stiffness properties for panels
        for i=1:oldPanelNum
            A=size(panel0{i});
            panelNodeNum=A(2);
            panelThick=panelThicknessMat(i);
            A=panel0{i};

            lsum=0;
            area=0;
            for j=1:panelNodeNum
                barNum1=SearchCreaseNum(barConnect,...
                    centerNodeStart+i,A(j));
                node1=newNode(centerNodeStart+i,:);
                node2=newNode(A(j),:);
                if j==panelNodeNum
                    barNum2=SearchCreaseNum(barConnect,...
                        A(1),A(j));
                    node3=newNode(A(1),:);
                else
                    barNum2=SearchCreaseNum(barConnect,...
                        A(j),A(j+1));
                    node3=newNode(A(j+1),:);                       
                end
                lsum=lsum+barLength(barNum1);
                lsum=lsum+barLength(barNum2);

                area=area+0.5*norm(cross(node1-node2,node3-node2));
            end
            aPanel=CalculateApanel(panelThick,area,panelPoisson,lsum);
            kSprPanel=CalculateKsprpanel(panelW,panelThick,panelE);

            
            for j=1:panelNodeNum
                barNum1=SearchCreaseNum(barConnect,...
                    centerNodeStart+i,A(j));
                if j==panelNodeNum
                    barNum2=SearchCreaseNum(barConnect,...
                        A(1),A(j));
                else
                    barNum2=SearchCreaseNum(barConnect,...
                        A(j),A(j+1));                      
                end
                barArea(barNum1)=barArea(barNum1)+aPanel;
                barArea(barNum2)=barArea(barNum1)+aPanel;
                sprK(barNum1)=...
                    kSprPanel*barLength(barNum1);
            end
        end        
    end
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