%% Addjust stiffness matrix K to consider support
%
% This code exclude the corresponding rows and columns in the support by
% putting zeros in the corresponding position.
% This gives the exact solution
%
% Input: 
%       K: stiffness matrix before considering support;
%       supp: information about support;
%       Tinput: forces before considering support;
%       elasticSupportOpen: indicate if non-rigid support is used;
%       suppElastic: non-rigid support information;
%       U: current diformation field;
% Output:
%       Kwsupp: modified stiffness matrix with support;
%       T: modified forces with support;
%

function [Kwsupp,T]=Solver_ModKforSupp(obj,K,supp,Tinput,elasticSupportOpen,suppElastic,U)
    Kwsupp=K;
    A=size(supp);
    SuppSize=A(1);
    A=size(K);
    N=A(1);
    T=Tinput;
    for i=1:SuppSize
        TempNodeNum=supp(i,1);
        if supp(i,2)==1
            Kvv=K(TempNodeNum*3-2,TempNodeNum*3-2);
            Kwsupp(TempNodeNum*3-2,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3-2)=zeros(N,1);
            Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=Kvv;
            T(TempNodeNum*3-2)=0;
        end
        if supp(i,3)==1
            Kvv=K(TempNodeNum*3-1,TempNodeNum*3-1);
            Kwsupp(TempNodeNum*3-1,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3-1)=zeros(N,1);
            Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=Kvv;
            T(TempNodeNum*3-1)=0;
        end
        if supp(i,4)==1
            Kvv=K(TempNodeNum*3,TempNodeNum*3);
            Kwsupp(TempNodeNum*3,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3)=zeros(N,1);
            Kwsupp(TempNodeNum*3,TempNodeNum*3)=Kvv;
            T(TempNodeNum*3)=0;
        end    
    end
    
    ElaSuppSize=size(suppElastic,1);
    if elasticSupportOpen==1
        for i=1:ElaSuppSize
            TempNodeNum=suppElastic(i,1);
            Kwsupp((TempNodeNum-1)*3+suppElastic(i,2),(TempNodeNum-1)*3+suppElastic(i,2))= ...
                Kwsupp((TempNodeNum-1)*3+suppElastic(i,2),(TempNodeNum-1)*3+suppElastic(i,2))+suppElastic(i,3);
            T((TempNodeNum-1)*3+suppElastic(i,2))=...
                T((TempNodeNum-1)*3+suppElastic(i,2))...
                -suppElastic(i,3)*U(TempNodeNum,suppElastic(i,2));
        end      
    end
end