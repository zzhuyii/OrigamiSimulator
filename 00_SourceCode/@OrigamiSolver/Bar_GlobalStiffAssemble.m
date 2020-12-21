%% Assemble stiffness of bars
% This function assembles the global stiffness matrix of the bars
%
% Input: 
%       U: current deformation field;
%       Sx: the second PK stress for bar elements;
%       C: tangent stiffness of bar elements;
%       barArea: area of bars;
%       barLength: length of bars;
%       barConnect: the two nodal number of bars;
%       newNode: nodal coordinates of meshed origami;
% Output: 
%       Kbar: global stiffness matrix from bars;
%

function [Kbar]=Bar_GlobalStiffAssemble(obj,U,Sx,C,barArea,barLength,barConnect,newNode)
    A=size(newNode);
    N=A(1);
    Kbar=zeros(3*N,3*N);
    A=size(C);
    N=A(1);
    for i=1:N
       NodeIndex1=barConnect(i,1);
       NodeIndex2=barConnect(i,2);
       node1=newNode(NodeIndex1,:);
       node2=newNode(NodeIndex2,:);

       B1=1/(barLength(i)^2)*[-(node2-node1) (node2-node1)];
       iden=eye(3);
       B2=1/(barLength(i)^2)*[iden -iden; -iden iden];

       Utemp=[U(NodeIndex1,:)';U(NodeIndex2,:)'];   
       Ktemp=C(i)*barArea(i)*barLength(i)*(B1'+B2*Utemp)*(B1'+B2*Utemp)'+ ...
           Sx(i)*barArea(i)*barLength(i)*B2;

       index1=3*NodeIndex1-2;
       index2=3*NodeIndex2-2;

       Kbar(index1:index1+2,index1:index1+2)=Kbar(index1:index1+2,index1:index1+2)+Ktemp(1:3,1:3);
       Kbar(index2:index2+2,index2:index2+2)=Kbar(index2:index2+2,index2:index2+2)+Ktemp(4:6,4:6);
       Kbar(index1:index1+2,index2:index2+2)=Kbar(index1:index1+2,index2:index2+2)+Ktemp(1:3,4:6);
       Kbar(index2:index2+2,index1:index1+2)=Kbar(index2:index2+2,index1:index1+2)+Ktemp(4:6,1:3);
    end
end