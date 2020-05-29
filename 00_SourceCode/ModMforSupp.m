%% Addjust stiffness matrix K to consider support
% This code exclude the corresponding rows and columns in the support by
% putting zeros in the corresponding position.
% This gives the exact solution

function [Mwsupp]=ModMforSupp(M,Supp,ElasticSupportOpen,SuppElastic,U)
    Mwsupp=M;
    A=size(Supp);
    SuppSize=A(1);
    A=size(M);
    N=A(1);
    for i=1:SuppSize
        TempNodeNum=Supp(i,1);
        if Supp(i,2)==1
            Mvv=0;
            Mwsupp(TempNodeNum*3-2,1:N)=zeros(1,N);
            Mwsupp(1:N,TempNodeNum*3-2)=zeros(N,1);
            Mwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=Mvv;

        end
        if Supp(i,3)==1
            Mvv=0;
            Mwsupp(TempNodeNum*3-1,1:N)=zeros(1,N);
            Mwsupp(1:N,TempNodeNum*3-1)=zeros(N,1);
            Mwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=Mvv;

        end
        if Supp(i,4)==1
            Mvv=0;
            Mwsupp(TempNodeNum*3,1:N)=zeros(1,N);
            Mwsupp(1:N,TempNodeNum*3)=zeros(N,1);
            Mwsupp(TempNodeNum*3,TempNodeNum*3)=Mvv;

        end    
    end
    
    ElaSuppSize=size(SuppElastic,1);
    if ElasticSupportOpen==1
        for i=1:ElaSuppSize
            TempNodeNum=SuppElastic(i,1);
            Mwsupp((TempNodeNum-1)*3+SuppElastic(i,2),(TempNodeNum-1)*3+SuppElastic(i,2))= ...
                Mwsupp((TempNodeNum-1)*3+SuppElastic(i,2),(TempNodeNum-1)*3+SuppElastic(i,2))+SuppElastic(i,3);
        end      
    end
end