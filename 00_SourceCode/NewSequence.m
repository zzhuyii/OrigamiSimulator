%% Rearrange Numbering System for Sparse Matrix
% This code rearrange the numbering system so that the stiffness matrix can
% become a sparse matrix. 

% inverseNumbering: used to assemble stiffness matrix
% newNumbering: used to convert results back to the original numbering so
% that plotting can be performed easily

function [newNumbering,inverseNumbering]=NewSequence(Node)

    A=size(Node);
    NodeNum=A(1);
    NodeDistance=zeros(NodeNum,1);
    newNumbering=zeros(NodeNum,1);
    for i=1:NodeNum
       NodeDistance(i)=norm(Node(i,:)-Node(1,:)); 
    end
    [A,newNumbering]=sort(NodeDistance);
    inverseNumbering=zeros(NodeNum,1);
    for i=1:NodeNum
       inverseNumbering(newNumbering(i))=i; 
    end
    

    
    % The following code disable the code
    for i=1:NodeNum
        inverseNumbering(i)=i;
        newNumbering(i)=i;
    end
end