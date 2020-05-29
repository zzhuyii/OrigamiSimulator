%% Assembling of thermal conductivity matrix
% created: Yi Zhu 2020-04-20

function [ThermalMat,NodeNum]=...
    ThermalConductAssembleMat(newNode,newPanel,BarConnect,...
    BarLength,BarType,ModelConstant,oldPanel,Utemp,ThermalBCpanels)

    A=size(newPanel);
    PanelNum=A(2);
    A=size(BarConnect);
    BarNum=A(1);
    A=size(newNode);
    NodeNum=A(1);
    
    kpanel=ModelConstant{20};
    kcrease=ModelConstant{29};    
    kenv=ModelConstant{21};    
    tenv=ModelConstant{22};
    tenvCrease=ModelConstant{30};
    tpanel=ModelConstant{4};
    tcrease=ModelConstant{5};
    
    alpha=ModelConstant{32};% Thermal dissipation angle
    
    airLayer=ModelConstant{31};
    
    ThermalMat=zeros((airLayer+1)*NodeNum);
    maxEnvThickness=1000*ones(PanelNum,1);
    % This vector stores the maximum air layer thickness
    % Solve for maximum allowed air thickness

    % find the set of thermal BC panels
    BCindex=zeros(1);
    count=1;
    for i=1:PanelNum
        A=size(ThermalBCpanels);
        BCnum=A(1);
        check=0;
        for j=1:BCnum
            if oldPanel(i)==ThermalBCpanels(j)
                check=1;
            end
        end
        if check==1
            BCindex(count)=i;
            count=count+1;
        end
    end
       
    range=setdiff(1:PanelNum,BCindex);
    for i=range
        nodeInfo=newPanel{i};
        nodei=nodeInfo(1);
        nodej=nodeInfo(2);
        nodek=nodeInfo(3);
        pt1=(newNode(nodei,:)+Utemp(nodei,:))';
        pt2=(newNode(nodej,:)+Utemp(nodej,:))';
        pt3=(newNode(nodek,:)+Utemp(nodek,:))';
        
        D1=100000;
        D2=100000;
        D3=100000;
        
        for j=BCindex
            pointInfo=newPanel{j};
            pointi=pointInfo(1);
            pointj=pointInfo(2);
            pointk=pointInfo(3);            
            point1=(newNode(pointi,:)+Utemp(pointi,:))';
            point2=(newNode(pointj,:)+Utemp(pointj,:))';
            point3=(newNode(pointk,:)+Utemp(pointk,:))';            
            [D1temp,zone,N1,N2]=P2TDistance(pt1,point1,point2,point3);
            if D1>D1temp
                D1=D1temp;
            end
        end
        
        for j=BCindex
            pointInfo=newPanel{j};
            pointi=pointInfo(1);
            pointj=pointInfo(2);
            pointk=pointInfo(3);            
            point1=(newNode(pointi,:)+Utemp(pointi,:))';
            point2=(newNode(pointj,:)+Utemp(pointj,:))';
            point3=(newNode(pointk,:)+Utemp(pointk,:))';            
            [D2temp,zone,N1,N2]=P2TDistance(pt2,point1,point2,point3);
            if D2>D2temp
                D2=D2temp;
            end
        end
        
        for j=BCindex
            pointInfo=newPanel{j};
            pointi=pointInfo(1);
            pointj=pointInfo(2);
            pointk=pointInfo(3);            
            point1=(newNode(pointi,:)+Utemp(pointi,:))';
            point2=(newNode(pointj,:)+Utemp(pointj,:))';
            point3=(newNode(pointk,:)+Utemp(pointk,:))';            
            [D3temp,zone,N1,N2]=P2TDistance(pt3,point1,point2,point3);
            if D3>D3temp
                D3=D3temp;
            end
        end
        
        maxEnvThickness(i)=mean([D1,D2,D3]);
        
    end
    
    % assemble the thermal conductivity matrix
    for i=1:PanelNum
       nodeInfo=newPanel{i};
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
       for j=1:BarNum
           if nodeS==BarConnect(j,1) && nodeB==BarConnect(j,2)
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
       for j=1:BarNum
           if nodeS==BarConnect(j,1) && nodeB==BarConnect(j,2)
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
       for j=1:BarNum
           if nodeS==BarConnect(j,1) && nodeB==BarConnect(j,2)
               bar3=j;
           end
       end
       
       % Calculate the system properties
       Lsum=BarLength(bar1)+BarLength(bar2)+BarLength(bar3);
       thick = tcrease;
       if BarType(bar1)==5 || BarType(bar2)==5 || BarType(bar3)==5
           thick=tpanel(oldPanel(i));
       end
       
       % T3 for assembling the thermal matrix
       point1=(newNode(nodei,:)+Utemp(nodei,:))';
       point2=(newNode(nodej,:)+Utemp(nodej,:))';
       point3=(newNode(nodek,:)+Utemp(nodek,:))';  
       vec1=point2-point1;
       vec2=point3-point1;
       % angle between bar1 and bar2
       angle=acos(dot(vec1,vec2)/norm(vec1)/norm(vec2));
       % relocate the nodes onto a 2D plane
       x1=0;y1=0;
       x2=BarLength(bar1);y2=0;
       x3=BarLength(bar2)*cos(angle);y3=BarLength(bar2)*sin(angle);
       At=((x1*y2-x2*y1)+(x2*y3-x3*y2)+(x3*y1-x1*y3))/2;
       % calculate constants for shape functions
       a1=(x2*y3-x3*y2);b1=(y2-y3);c1=(x3-x2);
       a2=(x3*y1-x1*y3);b2=(y3-y1);c2=(x1-x3);
       a3=(x1*y2-x2*y1);b3=(y1-y2);c3=(x2-x1);
       
       kmat=kcrease;
       if oldPanel(i)~=0
           kmat=kpanel(oldPanel(i));
       end
       tempK=kmat*thick/4/(At)*[b1^2+c1^2  b1*b2+c1*c2  b1*b3+c1*c3;
                          b1*b2+c1*c2   b2^2+c2^2    b2*b3+c2*c3;
                          b1*b3+c1*c3    b2*b3+c2*c3   b3^2+c3^2];
       indexarray=[nodei,nodej,nodek];
       ThermalMat(indexarray,indexarray)=ThermalMat(indexarray,indexarray)+tempK;  
       
       % original bar area calculation. This method preserves the total
       % volume but it turns out that it is doing badly when tackling skew
       % triangle panels. Thus, we switch to CST elements.
       
%        V=At*thick;
%        Across=V/Lsum;       
%        Lsquaresum=BarLength(bar1)^2+BarLength(bar2)^2+BarLength(bar3)^2;
%        Across1=V*BarLength(bar1)/Lsquaresum;
%        Across2=V*BarLength(bar2)/Lsquaresum;
%        Across3=V*BarLength(bar3)/Lsquaresum;
%        if oldPanel(i)==0
%            kij=kcrease*Across1/BarLength(bar1);
%            kik=kcrease*Across2/BarLength(bar2);
%            kkj=kcrease*Across3/BarLength(bar3);
%        else           
%            kij=kpanel(oldPanel(i))*Across1/BarLength(bar1);
%            kik=kpanel(oldPanel(i))*Across2/BarLength(bar2);
%            kkj=kpanel(oldPanel(i))*Across3/BarLength(bar3);
%        end
%        % assemble the global matrix
%        % inside the structure
%        ThermalMat(nodei,nodei)=ThermalMat(nodei,nodei)+kij+kik;
%        ThermalMat(nodej,nodej)=ThermalMat(nodej,nodej)+kij+kkj;
%        ThermalMat(nodek,nodek)=ThermalMat(nodek,nodek)+kik+kkj;       
%        ThermalMat(nodei,nodej)=ThermalMat(nodei,nodej)-kij;
%        ThermalMat(nodej,nodei)=ThermalMat(nodej,nodei)-kij;
%        ThermalMat(nodej,nodek)=ThermalMat(nodej,nodek)-kkj;
%        ThermalMat(nodek,nodej)=ThermalMat(nodek,nodej)-kkj;
%        ThermalMat(nodei,nodek)=ThermalMat(nodei,nodek)-kik;
%        ThermalMat(nodek,nodei)=ThermalMat(nodek,nodei)-kik;
       
       % assemble the global matrix for
       % connection to the environment
       for q=1:airLayer
           airThick=0;
           if oldPanel(i)==0
               if tenvCrease<maxEnvThickness(i)
                   airThick=tenvCrease;
               else
                   airThick=maxEnvThickness(i);
               end
           else
               if tenv(oldPanel(i))<maxEnvThickness(i)
                   airThick=tenv(oldPanel(i));
               else
                   airThick=maxEnvThickness(i);
               end
           end
           
           Ac=2/3*At;
           Ac=Ac+Lsum*(q-0.5)/airLayer*airThick*2/3*tan(alpha);
           kic=kenv*Ac*airLayer/airThick;           
           kjc=kic;
           kkc=kic;

           ThermalMat(q*NodeNum+nodei,q*NodeNum+nodei)=...
               ThermalMat(q*NodeNum+nodei,q*NodeNum+nodei)+kic;
           ThermalMat(q*NodeNum+nodej,q*NodeNum+nodej)=...
               ThermalMat(q*NodeNum+nodej,q*NodeNum+nodej)+kjc;
           ThermalMat(q*NodeNum+nodek,q*NodeNum+nodek)=...
               ThermalMat(q*NodeNum+nodek,q*NodeNum+nodek)+kkc;

           ThermalMat((q-1)*NodeNum+nodei,(q-1)*NodeNum+nodei)=...
               ThermalMat((q-1)*NodeNum+nodei,(q-1)*NodeNum+nodei)+kic;
           ThermalMat((q-1)*NodeNum+nodej,(q-1)*NodeNum+nodej)=...
               ThermalMat((q-1)*NodeNum+nodej,(q-1)*NodeNum+nodej)+kjc;
           ThermalMat((q-1)*NodeNum+nodek,(q-1)*NodeNum+nodek)=...
               ThermalMat((q-1)*NodeNum+nodek,(q-1)*NodeNum+nodek)+kkc;

           ThermalMat(q*NodeNum+nodei,(q-1)*NodeNum+nodei)=...
               ThermalMat(q*NodeNum+nodei,(q-1)*NodeNum+nodei)-kic;
           ThermalMat((q-1)*NodeNum+nodei,q*NodeNum+nodei)=...
               ThermalMat((q-1)*NodeNum+nodei,q*NodeNum+nodei)-kic;

           ThermalMat(q*NodeNum+nodej,(q-1)*NodeNum+nodej)=...
               ThermalMat(q*NodeNum+nodej,(q-1)*NodeNum+nodej)-kjc;
           ThermalMat((q-1)*NodeNum+nodej,q*NodeNum+nodej)=...
               ThermalMat((q-1)*NodeNum+nodej,q*NodeNum+nodej)-kjc;

           ThermalMat(q*NodeNum+nodek,(q-1)*NodeNum+nodek)=...
               ThermalMat(q*NodeNum+nodek,(q-1)*NodeNum+nodek)-kkc;
           ThermalMat((q-1)*NodeNum+nodek,q*NodeNum+nodek)=...
               ThermalMat((q-1)*NodeNum+nodek,q*NodeNum+nodek)-kkc;
       end
    end

end