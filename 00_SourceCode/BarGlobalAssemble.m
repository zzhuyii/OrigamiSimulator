%% Assemble stiffness of bars
% This function assembles the global stiffness matrix of the bars

% The content in this code is based on the open-access deformable orgami 
% simulator developed by K. Liu and G. H. Paulino 
% [1] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  

function [Kbar]=BarGlobalAssemble(U,Sx,C,BarArea,BarLength,BarConnect,newNode)
    A=size(newNode);
    N=A(1);
    Kbar=zeros(3*N,3*N);
    A=size(C);
    N=A(1);
    for i=1:N
       NodeIndex1=BarConnect(i,1);
       NodeIndex2=BarConnect(i,2);
       node1=newNode(NodeIndex1,:);
       node2=newNode(NodeIndex2,:);

       B1=1/(BarLength(i)^2)*[-(node2-node1) (node2-node1)];
       iden=eye(3);
       B2=1/(BarLength(i)^2)*[iden -iden; -iden iden];

       Utemp=[U(NodeIndex1,:)';U(NodeIndex2,:)'];   
       Ktemp=C(i)*BarArea(i)*BarLength(i)*(B1'+B2*Utemp)*(B1'+B2*Utemp)'+ ...
           Sx(i)*BarArea(i)*BarLength(i)*B2;

       index1=3*NodeIndex1-2;
       index2=3*NodeIndex2-2;

       Kbar(index1:index1+2,index1:index1+2)=Kbar(index1:index1+2,index1:index1+2)+Ktemp(1:3,1:3);
       Kbar(index2:index2+2,index2:index2+2)=Kbar(index2:index2+2,index2:index2+2)+Ktemp(4:6,4:6);
       Kbar(index1:index1+2,index2:index2+2)=Kbar(index1:index1+2,index2:index2+2)+Ktemp(1:3,4:6);
       Kbar(index2:index2+2,index1:index1+2)=Kbar(index2:index2+2,index1:index1+2)+Ktemp(4:6,1:3);
    end
end