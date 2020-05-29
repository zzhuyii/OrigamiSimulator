%% Addjust stiffness matrix K to consider support
% This code exclude the corresponding rows and columns in the support by
% putting zeros in the corresponding position.
% This gives the exact solution

function [Kwsupp,T]=ModKforSupp(K,Supp,Tinput,ElasticSupportOpen,SuppElastic,U)
    Kwsupp=K;
    A=size(Supp);
    SuppSize=A(1);
    A=size(K);
    N=A(1);
    T=Tinput;
    for i=1:SuppSize
        TempNodeNum=Supp(i,1);
        if Supp(i,2)==1
            Kvv=K(TempNodeNum*3-2,TempNodeNum*3-2);
            Kwsupp(TempNodeNum*3-2,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3-2)=zeros(N,1);
            Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=Kvv;
            T(TempNodeNum*3-2)=0;
        end
        if Supp(i,3)==1
            Kvv=K(TempNodeNum*3-1,TempNodeNum*3-1);
            Kwsupp(TempNodeNum*3-1,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3-1)=zeros(N,1);
            Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=Kvv;
            T(TempNodeNum*3-1)=0;
        end
        if Supp(i,4)==1
            Kvv=K(TempNodeNum*3,TempNodeNum*3);
            Kwsupp(TempNodeNum*3,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3)=zeros(N,1);
            Kwsupp(TempNodeNum*3,TempNodeNum*3)=Kvv;
            T(TempNodeNum*3)=0;
        end    
    end
    
    ElaSuppSize=size(SuppElastic,1);
    if ElasticSupportOpen==1
        for i=1:ElaSuppSize
            TempNodeNum=SuppElastic(i,1);
            Kwsupp((TempNodeNum-1)*3+SuppElastic(i,2),(TempNodeNum-1)*3+SuppElastic(i,2))= ...
                Kwsupp((TempNodeNum-1)*3+SuppElastic(i,2),(TempNodeNum-1)*3+SuppElastic(i,2))+SuppElastic(i,3);
            T((TempNodeNum-1)*3+SuppElastic(i,2))=...
                T((TempNodeNum-1)*3+SuppElastic(i,2))...
                -SuppElastic(i,3)*U(TempNodeNum,SuppElastic(i,2));
        end      
    end
end