%% Calculate rotation of springs
% the function calculate the rotation of springs when given the deformed
% configuration of the structure

% The content in this code is based on the open-access deformable orgami 
% simulator developed by K. Liu and G. H. Paulino 
% [1] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  

function [Theta]=CreaseTheta(U,CreaseIJKL,newNode)
    A=size(CreaseIJKL);
    N=A(1);
    Theta=zeros(N,1);
    for i=1:N
        if CreaseIJKL(i,1)==0
        else        
            nodei=newNode(CreaseIJKL(i,1),:)+U(CreaseIJKL(i,1),:);
            nodej=newNode(CreaseIJKL(i,2),:)+U(CreaseIJKL(i,2),:);
            nodek=newNode(CreaseIJKL(i,3),:)+U(CreaseIJKL(i,3),:);
            nodel=newNode(CreaseIJKL(i,4),:)+U(CreaseIJKL(i,4),:);
            rij=nodei-nodej;
            rkj=nodek-nodej;
            rkl=nodek-nodel;
            m=cross(rij,rkj);
            n=cross(rkj,rkl);
            yita=1;
            if dot(m,rkl)==0
                yita=1;
            else
                yita=sign(dot(m,rkl));
            end
            Theta(i)=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi);
        end
    end
end