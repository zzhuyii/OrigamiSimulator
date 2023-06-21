%% Connector Assemble Stiffness matrix
% This function calcualte the stiffness matrix of connector elements

function connectorKmat = Connector_Stiffness(obj,currentU,newNode,connectorK,connectorNode)

    A=size(connectorNode);
    connectorNum=A(1);
    
    A=size(currentU);
    NodeNum=A(1);
    
    connectorKmat=zeros(3*NodeNum,3*NodeNum);
    
    for i=1:connectorNum
        index1=connectorNode(i,1);
        index2=connectorNode(i,2);
        
        node1=newNode(index1,:)+currentU(index1,:);
        node2=newNode(index2,:)+currentU(index2,:);        

        connectorKmat(3*index1-2:3*index1,3*index1-2:3*index1)=connectorK*eye(3);
        connectorKmat(3*index2-2:3*index2,3*index1-2:3*index1)=-connectorK*eye(3);
        connectorKmat(3*index1-2:3*index1,3*index2-2:3*index2)=-connectorK*eye(3);
        connectorKmat(3*index2-2:3*index2,3*index2-2:3*index2)=connectorK*eye(3);
        
    end
end

