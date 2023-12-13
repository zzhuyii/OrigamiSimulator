%% Assembling of thermal capacity matrix

function [thermalCapacityMat]=Thermal_AssembleCapacityMat(obj,thermal,U)

    A=size(obj.newPanel);
    newPanelNum=A(2);
    A=size(obj.barConnect);
    barNum=A(1);
    A=size(obj.newNode);
    thermalNodeNum=A(1);
   
    thermalCapacityMat=zeros((obj.envLayer+1)*thermalNodeNum);
    maxEnvThickness=10000000000*ones(newPanelNum,1);
    % This vector stores the maximum air layer thickness
    % Solve for maximum allowed air thickness

    % find the set of thermal BC panels
    BCindex=[];
    count=1;
    
    for i=1:newPanelNum
        A=size(thermal.thermalBoundaryPanelVec);
        BCnum=A(1);
        check=0;
        for j=1:BCnum
            if obj.newPanel2OldPanel(i)==thermal.thermalBoundaryPanelVec(j)
                check=1;
            end
        end
        if check==1
            BCindex(count)=i;
            count=count+1;
        end
    end
       
    % identify the set of penals that are not thermal boundary
    range=setdiff(1:newPanelNum,BCindex);
    for i=range
        nodeInfo=obj.newPanel{i};
        nodei=nodeInfo(1);
        nodej=nodeInfo(2);
        nodek=nodeInfo(3);
        pt1=(obj.newNode(nodei,:)+U(nodei,:))';
        pt2=(obj.newNode(nodej,:)+U(nodej,:))';
        pt3=(obj.newNode(nodek,:)+U(nodek,:))';
        
        D1=100000;
        D2=100000;
        D3=100000;
        
        for j=BCindex
            pointInfo=obj.newPanel{j};
            pointi=pointInfo(1);
            pointj=pointInfo(2);
            pointk=pointInfo(3);            
            point1=(obj.newNode(pointi,:)+U(pointi,:))';
            point2=(obj.newNode(pointj,:)+U(pointj,:))';
            point3=(obj.newNode(pointk,:)+U(pointk,:))';            
            [D1temp,zone,N1,N2]=obj.Contact_P2TDistance(pt1,point1,point2,point3);
            if D1>D1temp
                D1=D1temp;
            end
        end
        
        for j=BCindex
            pointInfo=obj.newPanel{j};
            pointi=pointInfo(1);
            pointj=pointInfo(2);
            pointk=pointInfo(3);            
            point1=(obj.newNode(pointi,:)+U(pointi,:))';
            point2=(obj.newNode(pointj,:)+U(pointj,:))';
            point3=(obj.newNode(pointk,:)+U(pointk,:))';            
            [D2temp,zone,N1,N2]=obj.Contact_P2TDistance(pt2,point1,point2,point3);
            if D2>D2temp
                D2=D2temp;
            end
        end
        
        for j=BCindex
            pointInfo=obj.newPanel{j};
            pointi=pointInfo(1);
            pointj=pointInfo(2);
            pointk=pointInfo(3);            
            point1=(obj.newNode(pointi,:)+U(pointi,:))';
            point2=(obj.newNode(pointj,:)+U(pointj,:))';
            point3=(obj.newNode(pointk,:)+U(pointk,:))';            
            [D3temp,zone,N1,N2]=obj.Contact_P2TDistance(pt3,point1,point2,point3);
            if D3>D3temp
                D3=D3temp;
            end
        end
        
        maxEnvThickness(i)=mean([D1,D2,D3]);
        
    end
    
    % assemble the thermal capacity matrix for triangle panels
    for i=1:newPanelNum
        
       nodeInfo=obj.newPanel{i};
       nodei=nodeInfo(1);
       nodej=nodeInfo(2);
       nodek=nodeInfo(3);
       
%        % identify bar numbers
%        % bar1 connects node i and node j
%        if nodei<nodej
%            nodeS=nodei;nodeB=nodej;
%        else
%            nodeS=nodej;nodeB=nodei;
%        end
%        bar1=0;
%        for j=1:barNum
%            if nodeS==obj.barConnect(j,1) && nodeB==obj.barConnect(j,2)
%                bar1=j;
%            end
%        end
%        % bar2 connects node i and node k
%        if nodei<nodek
%            nodeS=nodei;nodeB=nodek;
%        else
%            nodeS=nodek;nodeB=nodei;
%        end       
%        bar2=0;
%        for j=1:barNum
%            if nodeS==obj.barConnect(j,1) && nodeB==obj.barConnect(j,2)
%                bar2=j;
%            end
%        end
%        % bar3 connects node j and node k
%        if nodek<nodej
%            nodeS=nodek;nodeB=nodej;
%        else
%            nodeS=nodej;nodeB=nodek;
%        end       
%        bar3=0;
%        for j=1:barNum
%            if nodeS==obj.barConnect(j,1) && nodeB==obj.barConnect(j,2)
%                bar3=j;
%            end
%        end
       
       bar1=obj.newPanel2NewBar(i,1);
       bar2=obj.newPanel2NewBar(i,2);
       bar3=obj.newPanel2NewBar(i,3);

       % Calculate the system properties
       Lsum=obj.barLength(bar1)+obj.barLength(bar2)+obj.barLength(bar3);
       thick=0;
       if obj.barType(bar1)==5 || obj.barType(bar2)==5 || obj.barType(bar3)==5
           thick=obj.panelThickVec(obj.newPanel2OldPanel(i));
       else
           thick = obj.creaseThickVec(max(max(obj.newCrease2OldCrease(bar1),...
               obj.newCrease2OldCrease(bar2)),obj.newCrease2OldCrease(bar3)));
       end
       
       % Find the three nodes of triangle panels
       point1=(obj.newNode(nodei,:)+U(nodei,:))';
       point2=(obj.newNode(nodej,:)+U(nodej,:))';
       point3=(obj.newNode(nodek,:)+U(nodek,:))';  
       vec1=point2-point1;
       vec2=point3-point1;

       % angle between bar1 and bar2
       angle=acos(dot(vec1,vec2)/norm(vec1)/norm(vec2));

       % relocate the nodes onto a 2D plane
       x1=0;y1=0;
       x2=obj.barLength(bar1);y2=0;
       x3=obj.barLength(bar2)*cos(angle);y3=obj.barLength(bar2)*sin(angle);

       % area of triangle
       At=((x1*y2-x2*y1)+(x2*y3-x3*y2)+(x3*y1-x1*y3))/2;

       
       Cmat=obj.creaseThermalCapacity;
       density=obj.densityCrease;
       if obj.newPanel2OldPanel(i)~=0
           Cmat=obj.panelThermalConductVec(obj.newPanel2OldPanel(i));
           density=obj.densityPanel;
       end
       

       tempC=Cmat*thick/3*(At)*eye(3)*density;

       indexarray=[nodei,nodej,nodek];
       thermalCapacityMat(indexarray,indexarray)=thermalCapacityMat(indexarray,indexarray)+tempC;  
       
       % assemble the global capacity matrix for the 1D
       % connection bar elements to the environment
       for q=1:obj.envLayer
           airThick=0;     

%            if obj.newPanel2OldPanel(i)==0
%                if obj.t2RT<maxEnvThickness(i)
%                    airThick=obj.t2RT;
%                else
%                    airThick=maxEnvThickness(i);
%                end
%            else
%                if obj.t2RT<maxEnvThickness(i)
%                    airThick=obj.t2RT;
%                else
%                    airThick=maxEnvThickness(i);
%                end
%            end

           if obj.t2RT<maxEnvThickness(i)
               airThick=obj.t2RT;
           else
               airThick=maxEnvThickness(i);
           end
           
           % Solve for the bar area
           % We lump volume of submerged environment from two sides of the 
           % panel to three bar elements representing the surrounding

           Ac=2/3*At;
           Ac=Ac+Lsum*(q-0.5)/obj.envLayer*airThick*2/3*tan(obj.thermalDissipation);
           Cic=obj.envThermalCapacity*obj.envDensity*Ac*airThick/obj.envLayer/2;           
           Cjc=Cic;
           Ckc=Cic;

           thermalCapacityMat(q*thermalNodeNum+nodei,q*thermalNodeNum+nodei)=...
               thermalCapacityMat(q*thermalNodeNum+nodei,q*thermalNodeNum+nodei)+Cic;
           thermalCapacityMat(q*thermalNodeNum+nodej,q*thermalNodeNum+nodej)=...
               thermalCapacityMat(q*thermalNodeNum+nodej,q*thermalNodeNum+nodej)+Cjc;
           thermalCapacityMat(q*thermalNodeNum+nodek,q*thermalNodeNum+nodek)=...
               thermalCapacityMat(q*thermalNodeNum+nodek,q*thermalNodeNum+nodek)+Ckc;

           thermalCapacityMat((q-1)*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)=...
               thermalCapacityMat((q-1)*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)+Cic;
           thermalCapacityMat((q-1)*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)=...
               thermalCapacityMat((q-1)*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)+Cjc;
           thermalCapacityMat((q-1)*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)=...
               thermalCapacityMat((q-1)*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)+Ckc;

%            thermalCapacityMat(q*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)=...
%                thermalCapacityMat(q*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)-Cic;
%            thermalCapacityMat((q-1)*thermalNodeNum+nodei,q*thermalNodeNum+nodei)=...
%                thermalCapacityMat((q-1)*thermalNodeNum+nodei,q*thermalNodeNum+nodei)-Cic;
% 
%            thermalCapacityMat(q*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)=...
%                thermalCapacityMat(q*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)-Cjc;
%            thermalCapacityMat((q-1)*thermalNodeNum+nodej,q*thermalNodeNum+nodej)=...
%                thermalCapacityMat((q-1)*thermalNodeNum+nodej,q*thermalNodeNum+nodej)-Cjc;
% 
%            thermalCapacityMat(q*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)=...
%                thermalCapacityMat(q*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)-Ckc;
%            thermalCapacityMat((q-1)*thermalNodeNum+nodek,q*thermalNodeNum+nodek)=...
%                thermalCapacityMat((q-1)*thermalNodeNum+nodek,q*thermalNodeNum+nodek)-Ckc;
       end
    end

end