%% This function calculates the Derivatives of D in Zone 1
% 
% Input:
%       Point: the nodal coordinates of the Point;
%       T1: the nearest nodal coordinates of the triangle;
%       d: the point to triangle distance d;
% Output:
%       Dd2x2: the Hessian Matrix of the point to triangle distance;
%       Ddx: the Gradient Vector of the point to traingle distance;
%

function [Dd2x2,Ddx]=Contact_DerivativeZone1(Point,T1,d)    
    Ddx=[1/d*(Point-T1);1/d*(T1-Point)];  
    
    iden33=[1 0 0;
            0 1 0;
            0 0 1;];
        
    Dd2x0x0=1/d*iden33-1/(d^3)*(Point-T1)*(Point-T1)';
    Dd2x1x1=1/d*iden33-1/(d^3)*(T1-Point)*(T1-Point)';
    Dd2x0x1=-1/d*iden33-1/(d^3)*(Point-T1)*(T1-Point)';
    Dd2x1x0=-1/d*iden33-1/(d^3)*(T1-Point)*(Point-T1)';
    
    Dd2x2=[Dd2x0x0 Dd2x0x1;
           Dd2x1x0 Dd2x1x1;];  
end