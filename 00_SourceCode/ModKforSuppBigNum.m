%% Addjust stiffness matrix K to consider support
% This code exclude the corresponding rows and columns in the support by
% putting a big number in the corresponding position. 
% This gives the approximated solution

function [Kwsupp]=ModKforSuppBigNum(K,Supp)
    Kwsupp=K;
    A=size(Supp);
    SuppSize=A(1);
    A=size(K);
    N=A(1);
    for i=1:SuppSize
        TempNodeNum=(Supp(i,1));
        if Supp(i,2)==1
            Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=(10^8)+Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2);
        end
        if Supp(i,3)==1
            Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=(10^8)+Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1);
        end
        if Supp(i,4)==1
            Kwsupp(TempNodeNum*3,TempNodeNum*3)=(10^8)+Kwsupp(TempNodeNum*3,TempNodeNum*3);
        end    
    end
end