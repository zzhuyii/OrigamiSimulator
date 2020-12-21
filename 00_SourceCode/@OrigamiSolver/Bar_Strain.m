%% Calculate the strain of bars
% The fuction calculate the strains of bars when given the deformation
% information of the structures.
%
% Input: 
%       U: current deformation field;
%       barArea: area of bars;
%       barLength: length of bars;
%       barConnect: the two nodal number of bars;
%       newNode: nodal coordinates of meshed origami;
% Output: 
%       Ex: the strain vector of the bar elements;
%

function [Ex]=Bar_Strain (obj,U,newNode,barArea,barConnect,barLength)
    Ex=zeros(size(barArea));
    A=size(barArea);
    N=A(1);
    for i=1:N
        node1=barConnect(i,1);
        node2=barConnect(i,2);
        iden=eye(3);
        tempU=[U(node1,:)';U(node2,:)'];

        B1=1/(barLength(i)^2)*[-(newNode(node2,:)-newNode(node1,:)) (newNode(node2,:)-newNode(node1,:)) ];
        B2=1/(barLength(i)^2)*[iden -iden;-iden iden];
        Ex(i)=B1*tempU+0.5*tempU'*B2*tempU;
    end
end