%% This function calculates the Derivatives of D in Zone 0
% 
% Input:
%       Point: the nodal coordinates of the Point;
%       T1,T2,T3: the three nodal coordinates of the triangle;
% Output:
%       Dd2x2: the Hessian Matrix of the point to triangle distance;
%       Ddx: the Gradient Vector of the point to traingle distance;
%

function [Dd2x2,Ddx]=Contact_DerivativeZone0(Point,T1,T2,T3)
    Dd2x2=zeros(12,12);
    Ddx=zeros(12,1);
    abc=cross((T2-T1),(T3-T1));
    
    V=dot(Point-T1,abc);
    A=dot(abc,abc);
    
    DAx0=[0 0 0]';
    DAx1=2*cross((T2-T3),abc);
    DAx2=2*cross((T3-T1),abc);
    DAx3=2*cross((T1-T2),abc);
    
    DAx=[DAx0;DAx1;DAx2;DAx3];
    
    zero33=zeros(3);
    iden=[1 0 0;
          0 1 0;
          0 0 1;];
      
    DA2x1x1=2*(dot(T2-T3,T2-T3)*iden+(T2-T3)*(T3-T2)');
    DA2x2x2=2*(dot(T3-T1,T3-T1)*iden+(T3-T1)*(T1-T3)');
    DA2x3x3=2*(dot(T1-T2,T1-T2)*iden+(T1-T2)*(T2-T1)');
    
    DA2x1x2=-2*(-dot(T3-T1,T2-T3)*iden+(T3-T1)*(T2-T3)'+skew(abc));
    DA2x2x3=-2*(-dot(T1-T2,T3-T1)*iden+(T1-T2)*(T3-T1)'+skew(abc));
    DA2x3x1=-2*(-dot(T2-T3,T1-T2)*iden+(T2-T3)*(T1-T2)'+skew(abc));

    DA2x1x3=-2*(-dot(T1-T2,T2-T3)*iden+(T1-T2)*(T2-T3)'-skew(abc));
    DA2x3x2=-2*(-dot(T3-T1,T1-T2)*iden+(T3-T1)*(T1-T2)'-skew(abc));
    DA2x2x1=-2*(-dot(T2-T3,T3-T1)*iden+(T2-T3)*(T3-T1)'-skew(abc));
    
    DA2x2=[zero33 zero33 zero33 zero33;
           zero33 DA2x1x1 DA2x1x2 DA2x1x3;
           zero33 DA2x2x1 DA2x2x2 DA2x2x3;
           zero33 DA2x3x1 DA2x3x2 DA2x3x3;];
       
    DVx0=abc;
    DVx1=cross(T2-T3,Point-T1)-abc;
    DVx2=cross(T3-T1,Point-T1);
    DVx3=cross(T1-T2,Point-T1);
    
    DVx=[DVx0;DVx1;DVx2;DVx3];
    
    DV2x2=[zero33 -skew(T2-T3) -skew(T3-T1) -skew(T1-T2);
           skew(T2-T3) zero33 -skew(T3-Point) -skew(T2-Point);
           skew(T3-T1) skew(T3-Point) zero33 -skew(T1-Point);
           skew(T1-T2) skew(T2-Point) skew(T1-Point) zero33];
       
    yita=sign(V);
       
    for i=1:12
        Ddx(i)=yita*(A^(-1/2))*DVx(i)-0.5*yita*V*(A^(-3/2))*DAx(i);
        for j=1:12            
            Dd2x2(i,j)=yita/sqrt(A)*DV2x2(i,j) ...
                -0.5*yita/sqrt((A^3))*DVx(i)*DAx(j)...
                -0.5*yita*V/sqrt(A^3)*DA2x2(i,j)...
                -0.5*yita/sqrt(A^3)*DAx(i)*DVx(j) ...
                +3/4*yita*V/sqrt(A^5)*DAx(i)*DAx(j); 
        end        
    end
end

%% Additonal Functions
function [A]=skew(v)
    A=[0 -v(3) v(2);
       v(3) 0 -v(1);
       -v(2) v(1) 0;];
end