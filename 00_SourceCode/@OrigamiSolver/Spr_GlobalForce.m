%% Assemble inner force vector for springs
% The output Tspr is the global inner force vector for rotational springs.
%
% Input: 
%       U: displacement field;
%       M: bending moment of spring element;
%       sprIJKL: stores the connectivity of each rotational springs;
%       sprKadj: spring stiffness after adjustment;
%       newNode: the nodal coordinates of each new node;
% Output:
%       Tspr: internal forces vector of springs;
%

function [Tspr]=Spr_GlobalForce(obj,U,M,sprIJKL,sprKadj,newNode)
    
    A=size(U);
    NodeNum=A(1);
    
    Tspr=zeros(3*NodeNum,1);
    
%     Tspr_1=zeros(3*NodeNum,1);   
%     

%     localT_his={};
%     tempNum=1;
    
    %% The non-vectorized version of the code
    
%     A=size(sprKadj);    
%     BarNum=A(1);
%     for i=1:BarNum
% 
%         if sprIJKL(i,1)==0   
%         else
%             nodei=newNode(sprIJKL(i,1),:)+U(sprIJKL(i,1),:);
%             nodej=newNode(sprIJKL(i,2),:)+U(sprIJKL(i,2),:);
%             nodek=newNode(sprIJKL(i,3),:)+U(sprIJKL(i,3),:);
%             nodel=newNode(sprIJKL(i,4),:)+U(sprIJKL(i,4),:);
% 
%             index1=3*sprIJKL(i,1)-2;
%             index2=3*sprIJKL(i,2)-2;
%             index3=3*sprIJKL(i,3)-2;
%             index4=3*sprIJKL(i,4)-2;
% 
%             rij=(nodei-nodej)';
%             rkj=(nodek-nodej)';
%             rkl=(nodek-nodel)';
%             m=cross(rij,rkj);
%             n=cross(rkj,rkl);
% 
%             parti=norm(rkj)/norm(m)/norm(m)*m;
%             partl=-norm(rkj)/norm(n)/norm(n)*n;
%             partj=((dot(rij,rkj)/norm(rkj)/norm(rkj))-1)*parti-dot(rkl,rkj)/norm(rkj)/norm(rkj)*partl;
%             partk=((dot(rkl,rkj)/norm(rkj)/norm(rkj))-1)*partl-dot(rij,rkj)/norm(rkj)/norm(rkj)*parti;
% 
%             part=[parti;partj;partk;partl];  
%             localT=M(i)*part;    
% 
%             Tspr(index1:index1+2)=Tspr(index1:index1+2)+localT(1:3);
%             Tspr(index2:index2+2)=Tspr(index2:index2+2)+localT(4:6);
%             Tspr(index3:index3+2)=Tspr(index3:index3+2)+localT(7:9);
%             Tspr(index4:index4+2)=Tspr(index4:index4+2)+localT(10:12);
% 
%         end
%     end
     
    %% The vectorized version of the code
    
    spr_i=sprIJKL(:,1);
    spr_j=sprIJKL(:,2);
    spr_k=sprIJKL(:,3);
    spr_l=sprIJKL(:,4);

    nonzero=find(spr_i);
    spr_Num=length(nonzero);

    spr_i=nonzeros(spr_i);
    spr_j=nonzeros(spr_j);
    spr_k=nonzeros(spr_k);
    spr_l=nonzeros(spr_l);

    nodei=newNode(spr_i,:)+U(spr_i,:);
    nodej=newNode(spr_j,:)+U(spr_j,:);
    nodek=newNode(spr_k,:)+U(spr_k,:);
    nodel=newNode(spr_l,:)+U(spr_l,:);

    rij=(nodei-nodej);
    rkj=(nodek-nodej);
    rkl=(nodek-nodel);

    m=cross(rij,rkj);
    n=cross(rkj,rkl);  

    m_square=dot(m',m')';
    n_square=dot(n',n')';

    rkj_square=dot(rkj',rkj')';
    rkj_norm=sqrt(rkj_square);

    parti=(rkj_norm./m_square).*m;
    partl=-rkj_norm./n_square.*n;
    partj=((dot(rij',rkj')'./rkj_norm./rkj_norm-1)).*parti-dot(rkl',rkj')'./rkj_norm./rkj_norm.*partl;
    partk=((dot(rkl',rkj')'./rkj_norm./rkj_norm-1)).*partl-dot(rij',rkj')'./rkj_norm./rkj_norm.*parti;

    gradient=[parti,partj,partk,partl];
    
    M_temp=M(nonzero);
    
    localT=M_temp.*gradient;        
    
    index1=3*spr_i-2;
    index2=3*spr_j-2;
    index3=3*spr_k-2;
    index4=3*spr_l-2;

    index=[index1, index1+1, index1+2,...
             index2, index2+1, index2+2,...
             index3, index3+1, index3+2,...
             index4, index4+1, index4+2];
         
    
%     t=10;
%     localT_his{t}
%     localT(t,:)
    
    index=index(:);
    localT=localT(:);
    
    for i=1:length(index)
        Tspr(index(i))=Tspr(index(i))+localT(i);        
    end
end
