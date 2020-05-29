%% Calculate point to triangle distance
% This code calculates the closest distance between the point and the
% triangle for simulating contact.

% The content in this code is based on the contact simulation proposed 
% in the PRSA paper:
% [1] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' PRSA. (submitted) 

function [D,zone,N1,N2]=P2TDistance(NPoint,Npt1,Npt2,Npt3)

    N1=0;
    N2=0;
    % Solve The Plane Equation
    a=(Npt2(2)-Npt1(2))*(Npt3(3)-Npt1(3))-(Npt2(3)-Npt1(3))*(Npt3(2)-Npt1(2));
    b=(Npt2(3)-Npt1(3))*(Npt3(1)-Npt1(1))-(Npt2(1)-Npt1(1))*(Npt3(3)-Npt1(3));
    c=(Npt2(1)-Npt1(1))*(Npt3(2)-Npt1(2))-(Npt2(2)-Npt1(2))*(Npt3(1)-Npt1(1));

    d=-Npt1(1)*(Npt2(2)-Npt1(2))*(Npt3(3)-Npt1(3)) ...
      -Npt1(2)*(Npt3(1)-Npt1(1))*(Npt2(3)-Npt1(3)) ...
      -Npt1(3)*(Npt2(1)-Npt1(1))*(Npt3(2)-Npt1(2)) ...
      +Npt1(3)*(Npt2(2)-Npt1(2))*(Npt3(1)-Npt1(1)) ...
      +Npt1(2)*(Npt2(1)-Npt1(1))*(Npt3(3)-Npt1(3)) ...
      +Npt1(1)*(Npt2(3)-Npt1(3))*(Npt3(2)-Npt1(2));

    P2PD=abs(a*NPoint(1)+b*NPoint(2)+c*NPoint(3)+d)/sqrt(a*a+b*b+c*c);
    planeNorm=[a;b;c];
    planeNorm=planeNorm/norm(planeNorm);

    B=Npt1;
    E0=Npt2-Npt1;
    E1=Npt3-Npt1;
    E2=cross(E0,E1);
    
    V=NPoint-B;
    K=[E0 E1 E2];
    Decomp=linsolve(K,V);  

   s0=Decomp(1);
    t0=Decomp(2);

    V1=cross(E2,E1);
    V2=cross(E0,E2);
    V3=cross(E1-E0,E2);

    if  s0>=0 && t0>=0 && s0+t0<=1 
        D=P2PD;
        zone=0;
    else
        Vec1=(NPoint-Npt1);
        K1=[V1 V2 cross(V1,V2)];
        Decomp=linsolve(K1,Vec1);

        s1=Decomp(1);
        t1=Decomp(2);

        Vec2=(NPoint-Npt2);
        K2=[V2 V3 cross(V2,V3)];
        Decomp=linsolve(K2,Vec2);

        s2=Decomp(1);
        t2=Decomp(2);

        Vec3=(NPoint-Npt3);
        K3=[V3 V1 cross(V3,V1)];
        Decomp=linsolve(K3,Vec3);

        s3=Decomp(1);
        t3=Decomp(2); 

        if s1>=0 && t1>=0
            D=norm(NPoint-Npt1);
            zone=1;
            N1=1;
        elseif s2>=0 && t2>=0
            D=norm(NPoint-Npt2);
            zone=1;
            N1=2;
        elseif s3>=0 && t3>=0
            D=norm(NPoint-Npt3);
            zone=1;
            N1=3;
        else            
            judge11=(NPoint-Npt2)'*V2;
            judge12=(NPoint-Npt1)'*E0;
            judge13=(NPoint-Npt2)'*(-E0);
            
            judge21=(NPoint-Npt3)'*V3;
            judge22=(NPoint-Npt2)'*(E1-E0);
            judge23=(NPoint-Npt3)'*(E0-E1);
            
            judge31=(NPoint-Npt1)'*V1;
            judge32=(NPoint-Npt1)'*E1;
            judge33=(NPoint-Npt3)'*(-E1);
            
            
            [D1]=PtoL(NPoint,Npt1,Npt2);
            [D2]=PtoL(NPoint,Npt2,Npt3);
            [D3]=PtoL(NPoint,Npt3,Npt1);        
            if (judge11>=0 && judge12>=0)&& judge13>=0
                D=D1;
                zone=2;
                N1=1;
                N2=2;
            elseif (judge21>=0 && judge22>=0)&& judge23>=0
                D=D2;
                zone=2;
                N1=2;
                N2=3;
            elseif (judge31>=0 && judge32>=0)&& judge33>=0
                D=D3;
                zone=2;
                N1=3;
                N2=1;
            else
                D=100;
                zone=100;
                fprintf('    Point to Panel Distance Caculation May Have Failed \n');
%             else
%                 D=0;
%                 zone=0;
            end
        end
    end
end

function [d,Direct]= PtoL(pt, v1, v2)
    a = v1 - v2;
    b = pt - v2;
    d = norm(cross(a,b)) / norm(a);
    normB=norm(b);
    c=sqrt(normB^2-d^2);

    aDirec=a/norm(a);
    Project=v2+aDirec*c;

    Direct=pt-Project;
end



