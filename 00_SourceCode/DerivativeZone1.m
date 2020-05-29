%% This function calculate the Derivatives of D in Zone 1
% This is a point to point case, where the point to triangle distance is
% just the distance between two nodes;

% The content in this code is based on the contact simulation proposed 
% in the PRSA paper:
% [1] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' PRSA. (submitted)  

function [Dd2x2,Ddx]=DerivativeZone1(Point,T1,d)    
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