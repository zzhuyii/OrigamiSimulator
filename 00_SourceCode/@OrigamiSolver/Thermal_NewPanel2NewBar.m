
%% This function calculate the NewPanel2NewBar matrix

function Thermal_NewPanel2NewBar(obj)
    
    A=size(obj.newPanel);
    newPanelNum=A(2);
    A=size(obj.barConnect);
    barNum=A(1);
    
    % assemble the thermal conductivity matrix
    obj.newPanel2NewBar=zeros(newPanelNum,3);
    
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
        
        obj.newPanel2NewBar(i,:)=[bar1,bar2,bar3];
    end
end

