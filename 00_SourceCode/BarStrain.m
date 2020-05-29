%% Calculate the strain of bars
% The fuction calculate the strains of bars when given the deformation
% information of the structures

% The content in this code is based on the open-access deformable orgami 
% simulator developed by K. Liu and G. H. Paulino 
% [1] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  

function [Ex]=BarStrain (U,Node,BarArea,BarConnect,BarLength)
    Ex=zeros(size(BarArea));
    A=size(BarArea);
    N=A(1);
    for i=1:N
        node1=BarConnect(i,1);
        node2=BarConnect(i,2);
        iden=eye(3);
        tempU=[U(node1,:)';U(node2,:)'];

        B1=1/(BarLength(i)^2)*[-(Node(node2,:)-Node(node1,:)) (Node(node2,:)-Node(node1,:)) ];
        B2=1/(BarLength(i)^2)*[iden -iden;-iden iden];
        Ex(i)=B1*tempU+0.5*tempU'*B2*tempU;
    end
end