%% Assemble global stiffness for springs
% This fucntion Calculate the global stiffness matrix from the 
% rotational spring elements
%
% Input: 
%       U: displacement field;
%       M: bending moment of spring element;
%       sprIJKL: stores the connectivity of each rotational springs;
%       sprKadj: spring stiffness after adjustment;
%       newNode: the nodal coordinates of each new node;
% Output:
%       Kspr: stiffness matrix of springs;
%

function [Kspr]=Spr_GlobalStiffAssemble(obj,U,M,sprIJKL,sprKadj,newNode)

    A=size(U);
    NodeNum=A(1);
    Kspr=zeros(3*NodeNum,3*NodeNum);
    
%     Kspr_1=zeros(3*NodeNum,3*NodeNum);
%     localK_record={};
%     part2_record={};
%     part2ii_record={};

%     tempNum=1;

%% This part of the code is not vectorized

%     A=size(sprKadj);
%     BarNum=A(1);
%     
%     for i=1:BarNum
% 
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
% 
%             m=cross(rij,rkj);
%             n=cross(rkj,rkl);     
% 
%             rkj_norm=norm(rkj);
%             m_square=m'*m;
%             n_square=n'*n;
% 
%             parti=norm(rkj)/m_square*m;
%             partl=-norm(rkj)/n_square*n;
% 
%             partj=((dot(rij,rkj)/rkj_norm/rkj_norm-1))*parti-dot(rkl,rkj)/rkj_norm/rkj_norm*partl;
%             partk=((dot(rkl,rkj)/rkj_norm/rkj_norm-1))*partl-dot(rij,rkj)/rkj_norm/rkj_norm*parti;
% 
%             part=[parti;partj;partk;partl];  
% 
%             A=dot(rij,rkj)/rkj_norm/rkj_norm;
%             B=dot(rkl,rkj)/rkj_norm/rkj_norm;
% 
%             partAj=1/rkj_norm/rkj_norm*((2*A-1)*rkj-rij);
%             partBj=1/rkj_norm/rkj_norm*(2*B*rkj-rkl);
%             partAk=1/rkj_norm/rkj_norm*(-2*A*rkj+rij);
%             partBk=1/rkj_norm/rkj_norm*((1-2*B)*rkj+rkl);
% 
%             part2ii=-norm(rkj)/(m_square*m_square)*(m*(cross(rkj,m))'+cross(rkj,m)*m');
%             part2ll=norm(rkj)/(n_square*n_square)*(n*(cross(rkj,n))'+cross(rkj,n)*n');
%             part2ik=m*rkj'/(m_square*norm(rkj))+norm(rkj)/(m_square*m_square)*(m*(cross(rij,m))'+cross(rij,m)*m');
%             part2lj=n*rkj'/(n_square*norm(rkj))-norm(rkj)/(n_square*n_square)*(n*(cross(rkl,n))'+cross(rkl,n)*n');
%             part2ij=-m*rkj'/(m_square*norm(rkj))+norm(rkj)/(m_square*m_square)*(m*(cross(rkj-rij,m))'+cross(rkj-rij,m)*m');
%             part2lk=-n*rkj'/(n_square*norm(rkj))-norm(rkj)/(n_square*n_square)*(n*(cross(rkj-rkl,n))'+cross(rkj-rkl,n)*n');
%             part2jj=parti*partAj'+(A-1)*part2ij-(partl*partBj'+B*part2lj);
%             part2jk=parti*partAk'+(A-1)*part2ik-(partl*partBk'+B*part2lk);
%             part2kk=partl*partBk'+(B-1)*part2lk-(parti*partAk'+A*part2ik);       
%             part2li=zeros(3);
% 
%             part2=[part2ii part2ij part2ik part2li';
%             part2ij' part2jj part2jk part2lj';
%             part2ik' part2jk' part2kk part2lk';
%             part2li part2lj part2lk part2ll];
%         
%             localK=sprKadj(i)*(part*part')+M(i)*part2;
% 
%             Kspr(index1:index1+2,index1:index1+2)=Kspr(index1:index1+2,index1:index1+2)+localK(1:3,1:3);
%             Kspr(index1:index1+2,index2:index2+2)=Kspr(index1:index1+2,index2:index2+2)+localK(1:3,4:6);
%             Kspr(index1:index1+2,index3:index3+2)=Kspr(index1:index1+2,index3:index3+2)+localK(1:3,7:9);
%             Kspr(index1:index1+2,index4:index4+2)=Kspr(index1:index1+2,index4:index4+2)+localK(1:3,10:12);
% 
%             Kspr(index2:index2+2,index1:index1+2)=Kspr(index2:index2+2,index1:index1+2)+localK(4:6,1:3);
%             Kspr(index2:index2+2,index2:index2+2)=Kspr(index2:index2+2,index2:index2+2)+localK(4:6,4:6);
%             Kspr(index2:index2+2,index3:index3+2)=Kspr(index2:index2+2,index3:index3+2)+localK(4:6,7:9);
%             Kspr(index2:index2+2,index4:index4+2)=Kspr(index2:index2+2,index4:index4+2)+localK(4:6,10:12);
% 
%             Kspr(index3:index3+2,index1:index1+2)=Kspr(index3:index3+2,index1:index1+2)+localK(7:9,1:3);
%             Kspr(index3:index3+2,index2:index2+2)=Kspr(index3:index3+2,index2:index2+2)+localK(7:9,4:6);
%             Kspr(index3:index3+2,index3:index3+2)=Kspr(index3:index3+2,index3:index3+2)+localK(7:9,7:9);
%             Kspr(index3:index3+2,index4:index4+2)=Kspr(index3:index3+2,index4:index4+2)+localK(7:9,10:12);
% 
%             Kspr(index4:index4+2,index1:index1+2)=Kspr(index4:index4+2,index1:index1+2)+localK(10:12,1:3);
%             Kspr(index4:index4+2,index2:index2+2)=Kspr(index4:index4+2,index2:index2+2)+localK(10:12,4:6);
%             Kspr(index4:index4+2,index3:index3+2)=Kspr(index4:index4+2,index3:index3+2)+localK(10:12,7:9);
%             Kspr(index4:index4+2,index4:index4+2)=Kspr(index4:index4+2,index4:index4+2)+localK(10:12,10:12);
% 
%         end 
%     end

%% This is the vectorized version of the code
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

    m=cross(rij,rkj,2);
    n=cross(rkj,rkl,2);  

    m_square=dot(m',m')';
    n_square=dot(n',n')';

    rkj_square=dot(rkj',rkj')';
    rkj_norm=sqrt(rkj_square);

    parti=(rkj_norm./m_square).*m;
    partl=-rkj_norm./n_square.*n;
    partj=((dot(rij',rkj')'./rkj_norm./rkj_norm-1)).*parti-dot(rkl',rkj')'./rkj_norm./rkj_norm.*partl;
    partk=((dot(rkl',rkj')'./rkj_norm./rkj_norm-1)).*partl-dot(rij',rkj')'./rkj_norm./rkj_norm.*parti;

    gradient=[parti,partj,partk,partl];

    A=dot(rij',rkj')'./rkj_norm./rkj_norm;
    B=dot(rkl',rkj')'./rkj_norm./rkj_norm;

    partAj=1./rkj_norm./rkj_norm.*((2.*A-1).*rkj-rij);
    partBj=1./rkj_norm./rkj_norm.*(2.*B.*rkj-rkl);
    partAk=1./rkj_norm./rkj_norm.*(-2.*A.*rkj+rij);
    partBk=1./rkj_norm./rkj_norm.*((1-2.*B).*rkj+rkl);   

    m1=m(:,1);
    m2=m(:,2);
    m3=m(:,3);
    n1=n(:,1);
    n2=n(:,2);
    n3=n(:,3);

    cross_rkj_m=cross(rkj',m')';
    part2ii_temp_1=m1.*cross(rkj',m')'+cross_rkj_m(:,1).*m;
    part2ii_temp_2=m2.*cross(rkj',m')'+cross_rkj_m(:,2).*m;
    part2ii_temp_3=m3.*cross(rkj',m')'+cross_rkj_m(:,3).*m;

    part2ii_temp=[part2ii_temp_1,part2ii_temp_2,part2ii_temp_3];
    part2ii_temp=-(rkj_norm./(m_square.*m_square)).*part2ii_temp;

    % part2ii=-rkj_norm./(m_square.*m_square).*(m*(cross(rkj',m')')+cross(rkj,m).*m');

    cross_rkj_n=cross(rkj',n')';
    part2ll_temp_1=n1.*cross(rkj',n')'+cross_rkj_n(:,1).*n;
    part2ll_temp_2=n2.*cross(rkj',n')'+cross_rkj_n(:,2).*n;
    part2ll_temp_3=n3.*cross(rkj',n')'+cross_rkj_n(:,3).*n;

    part2ll_temp=[part2ll_temp_1,part2ll_temp_2,part2ll_temp_3];
    part2ll_temp=(rkj_norm./(n_square.*n_square)).*part2ll_temp;

    % part2ll=rkj_norm./(n_square.*n_square).*(n*(cross(rkj',n')')+cross(rkj,n)*n');

    cross_rij_m=cross(rij',m')';
    part2ik_temp_1=m1.*rkj./m_square./rkj_norm+(rkj_norm./m_square./m_square).*((m1.*cross_rij_m)+(cross_rij_m(:,1).*m));
    part2ik_temp_2=m2.*rkj./m_square./rkj_norm+(rkj_norm./m_square./m_square).*((m2.*cross_rij_m)+(cross_rij_m(:,2).*m));
    part2ik_temp_3=m3.*rkj./m_square./rkj_norm+(rkj_norm./m_square./m_square).*((m3.*cross_rij_m)+(cross_rij_m(:,3).*m));

    part2ik_temp=[part2ik_temp_1,part2ik_temp_2,part2ik_temp_3];

    % part2ik=m*rkj'/(m_square*norm(rkj))+norm(rkj)/(m_square*m_square)*(m*(cross(rij,m))'+cross(rij,m)*m');

    cross_rkl_n=cross(rkl',n')';
    part2lj_temp_1=(n1.*rkj./n_square./rkj_norm)-(rkj_norm./n_square./n_square).*((n1.*cross_rkl_n)+(cross_rkl_n(:,1).*n));
    part2lj_temp_2=(n2.*rkj./n_square./rkj_norm)-(rkj_norm./n_square./n_square).*((n2.*cross_rkl_n)+(cross_rkl_n(:,2).*n));
    part2lj_temp_3=(n3.*rkj./n_square./rkj_norm)-(rkj_norm./n_square./n_square).*((n3.*cross_rkl_n)+(cross_rkl_n(:,3).*n));

    part2lj_temp=[part2lj_temp_1,part2lj_temp_2,part2lj_temp_3];

    % part2lj=n*rkj'/(n_square*norm(rkj))-norm(rkj)/(n_square*n_square)*(n*(cross(rkl,n))'+cross(rkl,n)*n');

    cross_rkj_rij_m=cross(rkj-rij,m);
    part2ij_temp_1=-m1.*rkj./(m_square.*rkj_norm)+rkj_norm./(m_square.*m_square).*((m1.*cross_rkj_rij_m)+(cross_rkj_rij_m(:,1).*m));
    part2ij_temp_2=-m2.*rkj./(m_square.*rkj_norm)+rkj_norm./(m_square.*m_square).*((m2.*cross_rkj_rij_m)+(cross_rkj_rij_m(:,2).*m));
    part2ij_temp_3=-m3.*rkj./(m_square.*rkj_norm)+rkj_norm./(m_square.*m_square).*((m3.*cross_rkj_rij_m)+(cross_rkj_rij_m(:,3).*m));

    part2ij_temp=[part2ij_temp_1,part2ij_temp_2,part2ij_temp_3];

    % part2ij=-m*rkj'/(m_square*norm(rkj))+norm(rkj)/(m_square*m_square)*(m*(cross(rkj-rij,m))'+cross(rkj-rij,m)*m');

    cross_rkj_rkl_n=cross(rkj-rkl,n);
    part2lk_temp_1=-(n1.*rkj./(n_square.*rkj_norm))-rkj_norm./(n_square.*n_square).*((n1.*cross_rkj_rkl_n)+(cross_rkj_rkl_n(:,1).*n));
    part2lk_temp_2=-(n2.*rkj./(n_square.*rkj_norm))-rkj_norm./(n_square.*n_square).*((n2.*cross_rkj_rkl_n)+(cross_rkj_rkl_n(:,2).*n));
    part2lk_temp_3=-(n3.*rkj./(n_square.*rkj_norm))-rkj_norm./(n_square.*n_square).*((n3.*cross_rkj_rkl_n)+(cross_rkj_rkl_n(:,3).*n));

    part2lk_temp=[part2lk_temp_1,part2lk_temp_2,part2lk_temp_3];

    %part2lk=-n*rkj'/(n_square*norm(rkj))-norm(rkj)/(n_square*n_square)*(n*(cross(rkj-rkl,n))'+cross(rkj-rkl,n)*n');

    parti_1=parti(:,1);
    parti_2=parti(:,2);
    parti_3=parti(:,3);

    partl_1=partl(:,1);
    partl_2=partl(:,2);
    partl_3=partl(:,3);

    part2jj_temp_1=(parti_1.*partAj)-(partl_1.*partBj);
    part2jj_temp_2=(parti_2.*partAj)-(partl_2.*partBj);
    part2jj_temp_3=(parti_3.*partAj)-(partl_3.*partBj);

    part2jj_temp=[part2jj_temp_1,part2jj_temp_2,part2jj_temp_3];
    part2jj_temp=part2jj_temp+(A-1).*part2ij_temp-B.*part2lj_temp;

    % part2jj=parti*partAj'+(A-1)*part2ij-(partl*partBj'+B*part2lj);

    part2jk_temp_1=(parti_1.*partAk)-(partl_1.*partBk);
    part2jk_temp_2=(parti_2.*partAk)-(partl_2.*partBk);
    part2jk_temp_3=(parti_3.*partAk)-(partl_3.*partBk);

    part2jk_temp=[part2jk_temp_1,part2jk_temp_2,part2jk_temp_3];
    part2jk_temp=part2jk_temp+(A-1).*part2ik_temp-B.*part2lk_temp;

    % part2jk=parti*partAk'+(A-1)*part2ik-(partl*partBk'+B*part2lk);

    part2kk_temp_1=(partl_1.*partBk)-(parti_1.*partAk);
    part2kk_temp_2=(partl_2.*partBk)-(parti_2.*partAk);
    part2kk_temp_3=(partl_3.*partBk)-(parti_3.*partAk);

    part2kk_temp=[part2kk_temp_1,part2kk_temp_2,part2kk_temp_3];
    part2kk_temp=part2kk_temp+(B-1).*part2lk_temp-A.*part2ik_temp;

    % part2kk=partl*partBk'+(B-1)*part2lk-(parti*partAk'+A*part2ik);

    part2li_temp=zeros(size(part2ii_temp));

    % part2li=zeros(3);

    part2ki_temp=part2ik_temp(:,[1,4,7,2,5,8,3,6,9]);
    part2jl_temp=part2lj_temp(:,[1,4,7,2,5,8,3,6,9]);
    part2ji_temp=part2ij_temp(:,[1,4,7,2,5,8,3,6,9]);
    part2kl_temp=part2lk_temp(:,[1,4,7,2,5,8,3,6,9]);
    part2kj_temp=part2jk_temp(:,[1,4,7,2,5,8,3,6,9]);
    part2il_temp=part2li_temp(:,[1,4,7,2,5,8,3,6,9]);

    sprKadj_temp=sprKadj(nonzero);
    M_temp=M(nonzero);

    localK_part1=[(sprKadj_temp.*(gradient(:,1)).*gradient),...
                  (sprKadj_temp.*(gradient(:,2)).*gradient),...
                  (sprKadj_temp.*(gradient(:,3)).*gradient),...
                  (sprKadj_temp.*(gradient(:,4)).*gradient),...
                  (sprKadj_temp.*(gradient(:,5)).*gradient),...
                  (sprKadj_temp.*(gradient(:,6)).*gradient),...
                  (sprKadj_temp.*(gradient(:,7)).*gradient),...
                  (sprKadj_temp.*(gradient(:,8)).*gradient),...
                  (sprKadj_temp.*(gradient(:,9)).*gradient),...
                  (sprKadj_temp.*(gradient(:,10)).*gradient),...
                  (sprKadj_temp.*(gradient(:,11)).*gradient),...
                  (sprKadj_temp.*(gradient(:,12)).*gradient),];
              
    part2ii_temp=reshape(part2ii_temp,[spr_Num,3,3]);
    part2ij_temp=reshape(part2ij_temp,[spr_Num,3,3]);
    part2ji_temp=reshape(part2ji_temp,[spr_Num,3,3]);
    part2ik_temp=reshape(part2ik_temp,[spr_Num,3,3]);
    part2ki_temp=reshape(part2ki_temp,[spr_Num,3,3]);
    part2il_temp=reshape(part2il_temp,[spr_Num,3,3]);
    part2li_temp=reshape(part2li_temp,[spr_Num,3,3]);
    
    part2jj_temp=reshape(part2jj_temp,[spr_Num,3,3]);
    part2jk_temp=reshape(part2jk_temp,[spr_Num,3,3]);
    part2kj_temp=reshape(part2kj_temp,[spr_Num,3,3]);
    part2jl_temp=reshape(part2jl_temp,[spr_Num,3,3]);
    part2lj_temp=reshape(part2lj_temp,[spr_Num,3,3]);
    
    part2kk_temp=reshape(part2kk_temp,[spr_Num,3,3]);
    part2kl_temp=reshape(part2kl_temp,[spr_Num,3,3]);
    part2lk_temp=reshape(part2lk_temp,[spr_Num,3,3]);
    
    part2ll_temp=reshape(part2ll_temp,[spr_Num,3,3]);
    
    hessian_1=cat(2,part2ii_temp,part2ij_temp,part2ik_temp,part2il_temp);
    hessian_2=cat(2,part2ji_temp,part2jj_temp,part2jk_temp,part2jl_temp);
    hessian_3=cat(2,part2ki_temp,part2kj_temp,part2kk_temp,part2kl_temp);
    hessian_4=cat(2,part2li_temp,part2lj_temp,part2lk_temp,part2ll_temp);
              
    hessian=cat(3,hessian_1,hessian_2,hessian_3,hessian_4);
         
%     hessian2=[part2ii_temp,part2ji_temp,part2ki_temp,part2li_temp,...
%              part2ij_temp,part2jj_temp,part2kj_temp,part2lj_temp,...
%              part2ik_temp,part2jk_temp,part2kk_temp,part2lk_temp,...
%              part2il_temp,part2jl_temp,part2kl_temp,part2ll_temp,];

    localK_part2=M_temp.*reshape(hessian,[spr_Num,144]);
            
    localK_n=localK_part1+localK_part2;
    % Finally, we get the stiffness matrix in the local coordinate


    % part2=[part2ii part2ij part2ik part2li';
    % part2ij' part2jj part2jk part2lj';
    % part2ik' part2jk' part2kk part2lk';
    % part2li part2lj part2lk part2ll];
    % 
    % localK=sprKadj(i)*(part*part')+M(i)*part2;

    
    % The final step is to make send the value of the local stiffness in to
    % the global stiffness using the index value;
    
    index1=3*spr_i-2;
    index2=3*spr_j-2;
    index3=3*spr_k-2;
    index4=3*spr_l-2;

    index_dim=[index1, index1+1, index1+2,...
             index2, index2+1, index2+2,...
             index3, index3+1, index3+2,...
             index4, index4+1, index4+2];
    
    index=zeros(spr_Num,12,12,2);
    
    index(:,1,:,2)=index_dim;    
    index(:,2,:,2)=index_dim;
    index(:,3,:,2)=index_dim;
    index(:,4,:,2)=index_dim;
    index(:,5,:,2)=index_dim;
    index(:,6,:,2)=index_dim;
    index(:,7,:,2)=index_dim;
    index(:,8,:,2)=index_dim;
    index(:,9,:,2)=index_dim;
    index(:,10,:,2)=index_dim;
    index(:,11,:,2)=index_dim;
    index(:,12,:,2)=index_dim;
    
    index(:,:,1,1)=index_dim;    
    index(:,:,2,1)=index_dim;
    index(:,:,3,1)=index_dim;
    index(:,:,4,1)=index_dim;
    index(:,:,5,1)=index_dim;
    index(:,:,6,1)=index_dim;
    index(:,:,7,1)=index_dim;
    index(:,:,8,1)=index_dim;
    index(:,:,9,1)=index_dim;
    index(:,:,10,1)=index_dim;
    index(:,:,11,1)=index_dim;
    index(:,:,12,1)=index_dim;
    
%     index=reshape(index, [spr_Num,144,2]);  
%     for i=1:spr_Num
%         for j=1:144
%             Kspr(index(i,j,1),index(i,j,2))=Kspr(index(i,j,1),index(i,j,2))+localK_n(i,j);
%         end
%     end
    
    index=reshape(index, [spr_Num*144,2]);
    localK_n=localK_n(:);
    
    for i=1:length(localK_n)
        Kspr(index(i,1),index(i,2))=Kspr(index(i,1),index(i,2))+localK_n(i);
    end

%     t=500
%     reshape(localK_n(t,:),[12,12])
%     localK_record{t}
%     
%     reshape(hessian(t,:),[12,12])
%     part2_record{t}
%     
%     reshape(part2ii_temp(t,:),[3,3])
%     part2ii_record{t}
    
end

