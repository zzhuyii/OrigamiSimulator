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
A=size(sprKadj);
BarNum=A(1);
for i=1:BarNum

    if sprIJKL(i,1)==0   
    else
        nodei=newNode(sprIJKL(i,1),:)+U(sprIJKL(i,1),:);
        nodej=newNode(sprIJKL(i,2),:)+U(sprIJKL(i,2),:);
        nodek=newNode(sprIJKL(i,3),:)+U(sprIJKL(i,3),:);
        nodel=newNode(sprIJKL(i,4),:)+U(sprIJKL(i,4),:);
        
        index1=3*sprIJKL(i,1)-2;
        index2=3*sprIJKL(i,2)-2;
        index3=3*sprIJKL(i,3)-2;
        index4=3*sprIJKL(i,4)-2;
        
        rij=(nodei-nodej)';
        rkj=(nodek-nodej)';
        rkl=(nodek-nodel)';
        m=cross(rij,rkj);
        n=cross(rkj,rkl);

        parti=norm(rkj)/norm(m)/norm(m)*m;
        partl=-norm(rkj)/norm(n)/norm(n)*n;
        partj=((dot(rij,rkj)/norm(rkj)/norm(rkj)-1))*parti-dot(rkl,rkj)/norm(rkj)/norm(rkj)*partl;
        partk=((dot(rkl,rkj)/norm(rkj)/norm(rkj)-1))*partl-dot(rij,rkj)/norm(rkj)/norm(rkj)*parti;

        part=[parti;partj;partk;partl];  
        
        A=dot(rij,rkj)/(norm(rkj)^2);
        B=dot(rkl,rkj)/(norm(rkj)^2);
        partAj=1/(norm(rkj)^2)*((2*A-1)*rkj-rij);
        partBj=1/(norm(rkj)^2)*(2*B*rkj-rkl);
        partAk=1/(norm(rkj)^2)*(-2*A*rkj+rij);
        partBk=1/(norm(rkj)^2)*((1-2*B)*rkj+rkl);        

        part2ii=-norm(rkj)/(norm(m)^4)*(m*(cross(rkj,m))'+cross(rkj,m)*m');
        part2ll=norm(rkj)/(norm(n)^4)*(n*(cross(rkj,n))'+cross(rkj,n)*n');
        part2ik=m*rkj'/((norm(m)^2)*norm(rkj))+norm(rkj)/(norm(m)^4)*(m*(cross(rij,m))'+cross(rij,m)*m');
        part2lj=n*rkj'/((norm(n)^2)*norm(rkj))-norm(rkj)/(norm(n)^4)*(n*(cross(rkl,n))'+cross(rkl,n)*n');
        part2ij=-m*rkj'/((norm(m)^2)*norm(rkj))+norm(rkj)/(norm(m)^4)*(m*(cross(rkj-rij,m))'+cross(rkj-rij,m)*m');
        part2lk=-n*rkj'/((norm(n)^2)*norm(rkj))-norm(rkj)/(norm(n)^4)*(n*(cross(rkj-rkl,n))'+cross(rkj-rkl,n)*n');
        part2jj=parti*partAj'+(A-1)*part2ij-(partl*partBj'+B*part2lj);
        part2jk=parti*partAk'+(A-1)*part2ik-(partl*partBk'+B*part2lk);
        part2kk=partl*partBk'+(B-1)*part2lk-(parti*partAk'+A*part2ik);        
        part2li=zeros(3);
        
        part2=[part2ii part2ij part2ik part2li';
               part2ij' part2jj part2jk part2lj';
               part2ik' part2jk' part2kk part2lk';
               part2li part2lj part2lk part2ll];
        
        localK=sprKadj(i)*(part*part')+M(i)*part2;
        
        Kspr(index1:index1+2,index1:index1+2)=Kspr(index1:index1+2,index1:index1+2)+localK(1:3,1:3);
        Kspr(index1:index1+2,index2:index2+2)=Kspr(index1:index1+2,index2:index2+2)+localK(1:3,4:6);
        Kspr(index1:index1+2,index3:index3+2)=Kspr(index1:index1+2,index3:index3+2)+localK(1:3,7:9);
        Kspr(index1:index1+2,index4:index4+2)=Kspr(index1:index1+2,index4:index4+2)+localK(1:3,10:12);
        
        Kspr(index2:index2+2,index1:index1+2)=Kspr(index2:index2+2,index1:index1+2)+localK(4:6,1:3);
        Kspr(index2:index2+2,index2:index2+2)=Kspr(index2:index2+2,index2:index2+2)+localK(4:6,4:6);
        Kspr(index2:index2+2,index3:index3+2)=Kspr(index2:index2+2,index3:index3+2)+localK(4:6,7:9);
        Kspr(index2:index2+2,index4:index4+2)=Kspr(index2:index2+2,index4:index4+2)+localK(4:6,10:12);
        
        Kspr(index3:index3+2,index1:index1+2)=Kspr(index3:index3+2,index1:index1+2)+localK(7:9,1:3);
        Kspr(index3:index3+2,index2:index2+2)=Kspr(index3:index3+2,index2:index2+2)+localK(7:9,4:6);
        Kspr(index3:index3+2,index3:index3+2)=Kspr(index3:index3+2,index3:index3+2)+localK(7:9,7:9);
        Kspr(index3:index3+2,index4:index4+2)=Kspr(index3:index3+2,index4:index4+2)+localK(7:9,10:12);
        
        Kspr(index4:index4+2,index1:index1+2)=Kspr(index4:index4+2,index1:index1+2)+localK(10:12,1:3);
        Kspr(index4:index4+2,index2:index2+2)=Kspr(index4:index4+2,index2:index2+2)+localK(10:12,4:6);
        Kspr(index4:index4+2,index3:index3+2)=Kspr(index4:index4+2,index3:index3+2)+localK(10:12,7:9);
        Kspr(index4:index4+2,index4:index4+2)=Kspr(index4:index4+2,index4:index4+2)+localK(10:12,10:12);
        
    end
end
