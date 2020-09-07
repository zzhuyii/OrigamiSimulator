%% Calculate point to triangle distance
%
% Input:
%       Point: the nodal coordinates of the Point;
%       T1,T2,T3: the three nodal coordinates of the triangle;
% Output:
%       d: the Point to Triangle distance;
%       zone: the zone the the point falls in;
%       N1: The node index of the closest node/ or one of the closest 
%           line segment;
%       N2: The node index of the seond node of the closest line segment;
%

function [d,zone,N1,N2]=Contact_P2TDistance(Point,T1,T2,T3)

    N1=0;
    N2=0;
    % Solve The Plane Equation
    a=(T2(2)-T1(2))*(T3(3)-T1(3))-(T2(3)-T1(3))*(T3(2)-T1(2));
    b=(T2(3)-T1(3))*(T3(1)-T1(1))-(T2(1)-T1(1))*(T3(3)-T1(3));
    c=(T2(1)-T1(1))*(T3(2)-T1(2))-(T2(2)-T1(2))*(T3(1)-T1(1));

    d=-T1(1)*(T2(2)-T1(2))*(T3(3)-T1(3)) ...
      -T1(2)*(T3(1)-T1(1))*(T2(3)-T1(3)) ...
      -T1(3)*(T2(1)-T1(1))*(T3(2)-T1(2)) ...
      +T1(3)*(T2(2)-T1(2))*(T3(1)-T1(1)) ...
      +T1(2)*(T2(1)-T1(1))*(T3(3)-T1(3)) ...
      +T1(1)*(T2(3)-T1(3))*(T3(2)-T1(2));

    P2PD=abs(a*Point(1)+b*Point(2)+c*Point(3)+d)/sqrt(a*a+b*b+c*c);
    planeNorm=[a;b;c];
    planeNorm=planeNorm/norm(planeNorm);

    B=T1;
    E0=T2-T1;
    E1=T3-T1;
    E2=cross(E0,E1);
    
    V=Point-B;
    K=[E0 E1 E2];
    Decomp=linsolve(K,V);  

   s0=Decomp(1);
    t0=Decomp(2);

    V1=cross(E2,E1);
    V2=cross(E0,E2);
    V3=cross(E1-E0,E2);

    if  s0>=0 && t0>=0 && s0+t0<=1 
        d=P2PD;
        zone=0;
    else
        Vec1=(Point-T1);
        K1=[V1 V2 cross(V1,V2)];
        Decomp=linsolve(K1,Vec1);

        s1=Decomp(1);
        t1=Decomp(2);

        Vec2=(Point-T2);
        K2=[V2 V3 cross(V2,V3)];
        Decomp=linsolve(K2,Vec2);

        s2=Decomp(1);
        t2=Decomp(2);

        Vec3=(Point-T3);
        K3=[V3 V1 cross(V3,V1)];
        Decomp=linsolve(K3,Vec3);

        s3=Decomp(1);
        t3=Decomp(2); 

        if s1>=0 && t1>=0
            d=norm(Point-T1);
            zone=1;
            N1=1;
        elseif s2>=0 && t2>=0
            d=norm(Point-T2);
            zone=1;
            N1=2;
        elseif s3>=0 && t3>=0
            d=norm(Point-T3);
            zone=1;
            N1=3;
        else            
            judge11=(Point-T2)'*V2;
            judge12=(Point-T1)'*E0;
            judge13=(Point-T2)'*(-E0);
            
            judge21=(Point-T3)'*V3;
            judge22=(Point-T2)'*(E1-E0);
            judge23=(Point-T3)'*(E0-E1);
            
            judge31=(Point-T1)'*V1;
            judge32=(Point-T1)'*E1;
            judge33=(Point-T3)'*(-E1);
            
            
            [D1]=PtoL(Point,T1,T2);
            [D2]=PtoL(Point,T2,T3);
            [D3]=PtoL(Point,T3,T1);        
            if (judge11>=0 && judge12>=0)&& judge13>=0
                d=D1;
                zone=2;
                N1=1;
                N2=2;
            elseif (judge21>=0 && judge22>=0)&& judge23>=0
                d=D2;
                zone=2;
                N1=2;
                N2=3;
            elseif (judge31>=0 && judge32>=0)&& judge33>=0
                d=D3;
                zone=2;
                N1=3;
                N2=1;
            else
                d=100;
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



