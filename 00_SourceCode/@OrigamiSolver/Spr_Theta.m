%% Calculate rotation of springs
% the function calculate the rotation of springs when given the deformed
% configuration of the structure
%
% Input: 
%       U: displacement field of the origami;
%       sprIJKL: the four node used to calculate the spring rotation;
%       newNode: the nodal coordinates of origami after meshing;
% Output:
%       theta: the rotaion angle of each spring elements;
%

function [theta]=Spr_Theta(obj,U,sprIJKL,newNode)
    A=size(sprIJKL);
    N=A(1);
    theta=zeros(N,1);
    for i=1:N
        if sprIJKL(i,1)==0
        else        
            nodei=newNode(sprIJKL(i,1),:)+U(sprIJKL(i,1),:);
            nodej=newNode(sprIJKL(i,2),:)+U(sprIJKL(i,2),:);
            nodek=newNode(sprIJKL(i,3),:)+U(sprIJKL(i,3),:);
            nodel=newNode(sprIJKL(i,4),:)+U(sprIJKL(i,4),:);
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
            theta(i)=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi);
        end
    end
end