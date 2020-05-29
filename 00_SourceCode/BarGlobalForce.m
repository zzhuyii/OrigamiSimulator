%% Assemble inner force of bars
% This function assembles the global inner force vector of bars 

% The content in this code is based on the open-access deformable orgami 
% simulator developed by K. Liu and G. H. Paulino 
% [1] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  

function [Tbar]=BarGlobalForce(U,Sx,C,BarArea,BarLength,BarConnect,newNode)
    A=size(newNode);
    N=A(1);
    Tbar=zeros(3*N,1);
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
       Ttemp=Sx(i)*BarArea(i)*BarLength(i)*(B1'+B2*Utemp);

       index1=3*NodeIndex1-2;
       index2=3*NodeIndex2-2;

       Tbar(index1:index1+2)=Tbar(index1:index1+2)+Ttemp(1:3);
       Tbar(index2:index2+2)=Tbar(index2:index2+2)+Ttemp(4:6);   
    end
end