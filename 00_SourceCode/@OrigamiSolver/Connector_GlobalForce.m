%% Connector Assemble Forces
% This function calcualte the internal forces of connector elements

function connectorF = Connector_GlobalForce(obj,currentU,newNode,connectorK,connectorNode)

    A=size(connectorNode);
    connectorNum=A(1);
    
    A=size(currentU);
    NodeNum=A(1);
    
    connectorF=zeros(3*NodeNum,1);
    
    for i=1:connectorNum
        index1=connectorNode(i,1);
        index2=connectorNode(i,2);
        
        node1=newNode(index1,:)+currentU(index1,:);
        node2=newNode(index2,:)+currentU(index2,:);        

        connectorF(3*index1-2:3*index1)=connectorK*(node1'-node2');
        connectorF(3*index2-2:3*index2)=connectorK*(node2'-node1');        
    end
end

