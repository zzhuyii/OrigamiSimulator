%% This function calculates the Derivatives of D for Zone 2
% 
% Input:
%       Point: the nodal coordinates of the Point;
%       T1,T2: the two nodal coordinates of the line segment that is  
%              closest to the poing on the triangle;
% Output:
%       Dd2x2: the Hessian Matrix of the point to triangle distance;
%       Ddx: the Gradient Vector of the point to traingle distance;
%

function [Dd2x2,Ddx]=Contact_DerivativeZone2(Point,T1,T2)
    Dd2x2=zeros(9,9);
    Ddx=zeros(9,1);
    
    abc=cross(T2-T1,Point-T1);
    a=abc(1);
    b=abc(2);
    c=abc(3);
    
    A=a^2+b^2+c^2;
    L=dot(T2-T1,T2-T1);
    
    DAx1=2*cross((T2-Point),abc);
    DAx2=2*cross((Point-T1),abc);
    DAx0=2*cross((T1-T2),abc);  
    
    DAx=[DAx0,DAx1,DAx2];
    
    iden33=[1 0 0;
            0 1 0;
            0 0 1;];        
        
    DA2x1x1=2*(dot(T2-Point,T2-Point)*iden33+(T2-Point)*(Point-T2)');
    DA2x2x2=2*(dot(Point-T1,Point-T1)*iden33+(Point-T1)*(T1-Point)');
    DA2x0x0=2*(dot(T1-T2,T1-T2)*iden33+(T1-T2)*(T2-T1)');
    
    DA2x1x2=-2*(-dot(Point-T1,T2-Point)*iden33+(Point-T1)*(T2-Point)'+skew(abc));
    DA2x2x0=-2*(-dot(T1-T2,Point-T1)*iden33+(T1-T2)*(Point-T1)'+skew(abc));
    DA2x0x1=-2*(-dot(T2-Point,T1-T2)*iden33+(T2-Point)*(T1-T2)'+skew(abc));

    DA2x1x0=-2*(-dot(T1-T2,T2-Point)*iden33+(T1-T2)*(T2-Point)'-skew(abc));
    DA2x0x2=-2*(-dot(Point-T1,T1-T2)*iden33+(Point-T1)*(T1-T2)'-skew(abc));
    DA2x2x1=-2*(-dot(T2-Point,Point-T1)*iden33+(T2-Point)*(Point-T1)'-skew(abc));
       
    DA2x2=[DA2x0x0 DA2x1x0 DA2x2x0;
           DA2x0x1 DA2x1x1 DA2x2x1;
           DA2x0x2 DA2x1x2 DA2x2x2;];
    
    DLx=[0;0;0; 2*(T1-T2);2*(T2-T1)];
    DL2x2=[zeros(3) zeros(3) zeros(3);
           zeros(3) 2*iden33 -2*iden33;
           zeros(3) -2*iden33 2*iden33;];
       
    for i=1:9
        Ddx(i)=0.5*1/sqrt(A*L)*DAx(i)-0.5*sqrt((A)/(L^3))*DLx(i);
        for j=1:9            
            Dd2x2(i,j)=0.5*1/sqrt(A*L)*DA2x2(i,j)-0.25*(1/sqrt((A^3)*L))*DAx(i)*DAx(j)...
                -0.25*1/sqrt((L^3)*A)*DAx(i)*DLx(j)-0.5*sqrt(A/(L^3))*DL2x2(i,j) ...
                -0.25*1/sqrt(A*(L^3))*DLx(i)*DAx(j)+0.75*sqrt(A)/sqrt(L^5)*DLx(i)*DLx(j);            
        end        
    end
    
end

%% Additional Functions
function [A]=skew(v)
    A=[0 -v(3) v(2);
       v(3) 0 -v(1);
       -v(2) v(1) 0;];
end