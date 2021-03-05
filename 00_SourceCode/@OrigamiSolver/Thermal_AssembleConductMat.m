%% Assembling of thermal conductivity matrix
% 
% this function assembles the thermal conductivity matrix for conducting
% analysis of the heat transfer problem. 
%
% Input:
%       newNode: nodal coordinates of the meshed origami;
%       newPanel: panel connectivity of the meshed origami;
%       barConnect: two nodes of each bar;
%       barLength: length of bar;
%       barType: the type of each bar element;
%       U: deformation field;
%       modelThermalConstant: stores heat tranfer inputs;
%       modelMechanicalConstant: stores mechanical inputs;
%       newCrease2OldCrease: mapping between new crease to old ones;
%       newPanel2OldPanel: mapping bewteen new panel to old panels;
% Output:
%       thermalMat: conducntivity matrix;
%       thermalNodeNum: total number of nodes (no environmental nodes);
%

function [thermalMat]=Thermal_AssembleConductMat(obj,thermal,U)

    A=size(obj.newPanel);
    newPanelNum=A(2);
    A=size(obj.barConnect);
    barNum=A(1);
    A=size(obj.newNode);
    thermalNodeNum=A(1);
   
    thermalMat=zeros((obj.envLayer+1)*thermalNodeNum);
    maxEnvThickness=10000000000*ones(newPanelNum,1);
    % This vector stores the maximum air layer thickness
    % Solve for maximum allowed air thickness

    % find the set of thermal BC panels
    BCindex=[];
    count=1;
    for i=1:newPanelNum
        A=size(thermal.thermalBoundaryPanelMat);
        BCnum=A(1);
        check=0;
        for j=1:BCnum
            if obj.newPanel2OldPanel(i)==thermal.thermalBoundaryPanelMat(j)
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
    
    % assemble the thermal conductivity matrix
    for i=1:newPanelNum
       nodeInfo=obj.newPanel{i};
       nodei=nodeInfo(1);
       nodej=nodeInfo(2);
       nodek=nodeInfo(3);
       
       % identify bar numbers
       % bar1 connects node i and node j
       if nodei<nodej
           nodeS=nodei;nodeB=nodej;
       else
           nodeS=nodej;nodeB=nodei;
       end
       bar1=0;
       for j=1:barNum
           if nodeS==obj.barConnect(j,1) && nodeB==obj.barConnect(j,2)
               bar1=j;
           end
       end
       % bar2 connects node i and node k
       if nodei<nodek
           nodeS=nodei;nodeB=nodek;
       else
           nodeS=nodek;nodeB=nodei;
       end       
       bar2=0;
       for j=1:barNum
           if nodeS==obj.barConnect(j,1) && nodeB==obj.barConnect(j,2)
               bar2=j;
           end
       end
       % bar3 connects node j and node k
       if nodek<nodej
           nodeS=nodek;nodeB=nodej;
       else
           nodeS=nodej;nodeB=nodek;
       end       
       bar3=0;
       for j=1:barNum
           if nodeS==obj.barConnect(j,1) && nodeB==obj.barConnect(j,2)
               bar3=j;
           end
       end
       
       % Calculate the system properties
       Lsum=obj.barLength(bar1)+obj.barLength(bar2)+obj.barLength(bar3);
       thick=0;
       if obj.barType(bar1)==5 || obj.barType(bar2)==5 || obj.barType(bar3)==5
           thick=obj.panelThickMat(obj.newPanel2OldPanel(i));
       else
           thick = obj.creaseThickMat(max(max(obj.newCrease2OldCrease(bar1),...
               obj.newCrease2OldCrease(bar2)),obj.newCrease2OldCrease(bar3)));
       end
       
       % T3 for assembling the thermal matrix
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
       At=((x1*y2-x2*y1)+(x2*y3-x3*y2)+(x3*y1-x1*y3))/2;
       % calculate constants for shape functions
       a1=(x2*y3-x3*y2);b1=(y2-y3);c1=(x3-x2);
       a2=(x3*y1-x1*y3);b2=(y3-y1);c2=(x1-x3);
       a3=(x1*y2-x2*y1);b3=(y1-y2);c3=(x2-x1);
       
       kmat=obj.creaseThermalConduct;
       if obj.newPanel2OldPanel(i)~=0
           kmat=obj.panelThermalConductMat(obj.newPanel2OldPanel(i));
       end
       tempK=kmat*thick/4/(At)*[b1^2+c1^2  b1*b2+c1*c2  b1*b3+c1*c3;
                          b1*b2+c1*c2   b2^2+c2^2    b2*b3+c2*c3;
                          b1*b3+c1*c3    b2*b3+c2*c3   b3^2+c3^2];
       indexarray=[nodei,nodej,nodek];
       thermalMat(indexarray,indexarray)=thermalMat(indexarray,indexarray)+tempK;  
       
       % assemble the global matrix for
       % connection to the environment
       for q=1:obj.envLayer
           airThick=0;
           if obj.newPanel2OldPanel(i)==0
               if obj.t2RT<maxEnvThickness(i)
                   airThick=obj.t2RT;
               else
                   airThick=maxEnvThickness(i);
               end
           else
               if obj.t2RT<maxEnvThickness(i)
                   airThick=obj.t2RT;
               else
                   airThick=maxEnvThickness(i);
               end
           end
           
           Ac=2/3*At;
           Ac=Ac+Lsum*(q-0.5)/obj.envLayer*airThick*2/3*tan(obj.thermalDissipation);
           kic=obj.envThermalConduct*Ac*obj.envLayer/airThick;           
           kjc=kic;
           kkc=kic;

           thermalMat(q*thermalNodeNum+nodei,q*thermalNodeNum+nodei)=...
               thermalMat(q*thermalNodeNum+nodei,q*thermalNodeNum+nodei)+kic;
           thermalMat(q*thermalNodeNum+nodej,q*thermalNodeNum+nodej)=...
               thermalMat(q*thermalNodeNum+nodej,q*thermalNodeNum+nodej)+kjc;
           thermalMat(q*thermalNodeNum+nodek,q*thermalNodeNum+nodek)=...
               thermalMat(q*thermalNodeNum+nodek,q*thermalNodeNum+nodek)+kkc;

           thermalMat((q-1)*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)=...
               thermalMat((q-1)*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)+kic;
           thermalMat((q-1)*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)=...
               thermalMat((q-1)*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)+kjc;
           thermalMat((q-1)*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)=...
               thermalMat((q-1)*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)+kkc;

           thermalMat(q*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)=...
               thermalMat(q*thermalNodeNum+nodei,(q-1)*thermalNodeNum+nodei)-kic;
           thermalMat((q-1)*thermalNodeNum+nodei,q*thermalNodeNum+nodei)=...
               thermalMat((q-1)*thermalNodeNum+nodei,q*thermalNodeNum+nodei)-kic;

           thermalMat(q*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)=...
               thermalMat(q*thermalNodeNum+nodej,(q-1)*thermalNodeNum+nodej)-kjc;
           thermalMat((q-1)*thermalNodeNum+nodej,q*thermalNodeNum+nodej)=...
               thermalMat((q-1)*thermalNodeNum+nodej,q*thermalNodeNum+nodej)-kjc;

           thermalMat(q*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)=...
               thermalMat(q*thermalNodeNum+nodek,(q-1)*thermalNodeNum+nodek)-kkc;
           thermalMat((q-1)*thermalNodeNum+nodek,q*thermalNodeNum+nodek)=...
               thermalMat((q-1)*thermalNodeNum+nodek,q*thermalNodeNum+nodek)-kkc;
       end
    end

end