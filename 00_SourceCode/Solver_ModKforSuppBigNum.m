%% Addjust stiffness matrix K to consider support
%
% This code exclude the corresponding rows and columns in the support by
% putting a big number in the corresponding position. 
% This gives the approximated solution, and I do not really use it in 
% common analysis. Sometime I use it to check the rank of stiffness matrix.
%
% Input: 
%       K: stiffness matrix before considering support;
%       supp: information about support;
% Output:
%       Kwsupp: modified stiffness matrix with support;
%


function [Kwsupp]=Solver_ModKforSuppBigNum(K,supp)
    Kwsupp=K;
    A=size(supp);
    SuppSize=A(1);
    A=size(K);
    N=A(1);
    for i=1:SuppSize
        TempNodeNum=(supp(i,1));
        if supp(i,2)==1
            Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=(10^8)+Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2);
        end
        if supp(i,3)==1
            Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=(10^8)+Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1);
        end
        if supp(i,4)==1
            Kwsupp(TempNodeNum*3,TempNodeNum*3)=(10^8)+Kwsupp(TempNodeNum*3,TempNodeNum*3);
        end    
    end
end