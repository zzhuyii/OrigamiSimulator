%%%%%%%%%%%%%%%%%%%%%%  Origami Contact Simulator  %%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Yi Zhu, and Evegueni T. Filipov
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% contact model presented in this simulation code. 
%
% Code Features: 
% (1) Simulate Contact in Origami patterns
% (2) Simulate Compliant Creases in Origami Patterns
% (3) Automated Meshing for Compliant Creases
% (4) Solver for both folding and loading available
% (5) Nonrigid support available
%
% Reference:
% [1] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' PRSA. (submitted)       
% [2] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [3] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  
% [4] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [5] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%%%%%%%%%%%%%%%%%%%  Origami Contact Simulator  %%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all;
close all;

N=100;
height=0.4;
step=zeros(N,1);
d=zeros(N,1);
f=zeros(N,1);
k=zeros(N,1);
maxdelta=1.0;

forceHis=zeros(N,4,3);
PointCord=zeros(N,3);

Dxx=zeros(N,1);
pt1his=zeros(N,3);
pt2his=zeros(N,3);
pt3his=zeros(N,3);

for i=1:N
    step(i)=-1+i*maxdelta/N;
    
    pt1=[0;0;0.2;];
    pt2=[1;0;0;];
    pt3=[0;1;0;];  
    
    pt1his(i,:)=pt1';
    pt2his(i,:)=pt2';
    pt3his(i,:)=pt3'; 
        
    Point=[-0.2+i*maxdelta/N;-0.2+i*maxdelta/N;height];
    limit=1;
    rate=1;
    
    tempF=zeros(12,1);
    tempK=zeros(12);
    
    localF1=zeros(6,1);
    localK1=zeros(6,6);
    
    localF2=zeros(9,1);
    localK2=zeros(9,9);

    [D,zone,N1,N2]=Contact_P2TDistance(Point,pt1,pt2,pt3);    
    d(i)=D;
    
    if zone==0
        0
        [Dd2x2,Ddx]=Contact_DerivativeZone0(Point,pt1,pt2,pt3);
        if D<limit
            phi=pi/2*(1-D/limit);
            phi0=-pi/2/limit;

            for i1=1:12
                tempF(i1)=tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1);
                for j1=1:12
                    tempK(i1,j1)=(tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                        -(phi0^2)*Ddx(i1)*Ddx(j1);
                end
            end
        end
    elseif zone==1
        1
        if N1==1
            pt=pt1;
        elseif N1==2
            pt=pt2;
        else
            pt=pt3;
        end
        [Dd2x2,Ddx]=Contact_DerivativeZone1(Point,pt,D);
        if D<limit
            phi=pi/2*(1-D/limit);
            phi0=-pi/2/limit;

            for i1=1:6
                localF1(i1)=tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1);
                for j1=1:6
                    localK1(i1,j1)=(tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                        -(phi0^2)*Ddx(i1)*Ddx(j1);
                end
            end
        end
        
        tempF(1:3)=localF1(1:3);
        tempK(1:3,1:3)=localK1(1:3,1:3);

        if N1==1
            tempF(4:6)=localF1(4:6);
            tempK(4:6,4:6)=localK1(4:6,4:6);
            tempK(1:3,4:6)=localK1(1:3,4:6);
            tempK(4:6,1:3)=localK1(4:6,1:3);
        elseif N1==2
            tempF(7:9)=localF1(4:6);
            tempK(7:9,7:9)=localK1(4:6,4:6);
            tempK(1:3,7:9)=localK1(1:3,4:6);
            tempK(7:9,1:3)=localK1(4:6,1:3);
        else
            tempF(10:12)=localF1(4:6);
            tempK(10:12,10:12)=localK1(4:6,4:6);
            tempK(1:3,10:12)=localK1(1:3,4:6);
            tempK(10:12,1:3)=localK1(4:6,1:3);
        end   
    else
        2
        if N1==1
            temppt1=pt1;
        elseif N1==2
            temppt1=pt2;
        else
            temppt1=pt3;
        end
        if N2==1
            temppt2=pt1;
        elseif N2==2
            temppt2=pt2;
        else
            temppt2=pt3;
        end
        [Dd2x2,Ddx]=Contact_DerivativeZone2(Point,temppt1,temppt2);
        if D<limit
            phi=pi/2*(1-D/limit);
            phi0=-pi/2/limit;

            for i1=1:9
                localF2(i1)=tan(phi)*phi0*Ddx(i1)-phi*phi0*Ddx(i1);
                for j1=1:9
                    localK2(i1,j1)=(tan(phi)-phi)*Dd2x2(i1,j1)*phi0+(phi0^2)*(sec(phi))^2*Ddx(i1)*Ddx(j1)...
                        -(phi0^2)*Ddx(i1)*Ddx(j1);
                end
            end
        end
        
        tempF(1:3)=localF2(1:3);
        tempK(1:3,1:3)=localK2(1:3,1:3);
        
        if N1==1
            tempF(4:6)=localF2(4:6);
            tempK(4:6,4:6)=localK2(4:6,4:6);
            tempK(1:3,4:6)=localK2(1:3,4:6);
            tempK(4:6,1:3)=localK2(4:6,1:3);
        elseif N1==2
            tempF(7:9)=localF2(4:6);
            tempK(7:9,7:9)=localK2(4:6,4:6);
            tempK(1:3,7:9)=localK2(1:3,4:6);
            tempK(7:9,1:3)=localK2(4:6,1:3);
        else
            tempF(10:12)=localF2(4:6);
            tempK(10:12,10:12)=localK2(4:6,4:6);
            tempK(1:3,10:12)=localK2(1:3,4:6);
            tempK(10:12,1:3)=localK2(4:6,1:3);
        end

        if N2==1
            tempF(4:6)=localF2(7:9);
            tempK(4:6,4:6)=localK2(7:9,7:9);
            tempK(1:3,4:6)=localK2(1:3,7:9);
            tempK(4:6,1:3)=localK2(7:9,1:3);
        elseif N2==2
            tempF(7:9)=localF2(7:9);
            tempK(7:9,7:9)=localK2(7:9,7:9);
            tempK(1:3,7:9)=localK2(1:3,7:9);
            tempK(7:9,1:3)=localK2(7:9,1:3);
        else
            tempF(10:12)=localF2(7:9);
            tempK(10:12,10:12)=localK2(7:9,7:9);
            tempK(1:3,10:12)=localK2(1:3,7:9);
            tempK(10:12,1:3)=localK2(7:9,1:3);
        end
    end  
   
    Dxx(i)=Dd2x2(5,5);
    Dxx(i)=Ddx(2);
    
    f(i)=tempF(2);
    k(i)=tempK(1,1);   
    
    forceHis(i,1,:)=tempF(1:3);
    forceHis(i,2,:)=tempF(4:6);
    forceHis(i,3,:)=tempF(7:9);
    forceHis(i,4,:)=tempF(10:12);
    
%     forceHis(i,1,:)=Ddx(1:3);
%     forceHis(i,2,:)=Ddx(4:6);
%     forceHis(i,3,:)=Ddx(7:9);
%     forceHis(i,4,:)=Ddx(10:12);
    
    PointCord(i,:)=Point(1:3)';
    
end

F=zeros(4,3);
for j=1:3
   F(j,:)=tempF(3*j+1:3*j+3)';
end

figure
plot(step,d);
figure
plot(step,f);
figure
plot(step,k);

figure
View1=45;
View2=15;
Vsize=1.5;

rate=0.3;

pauseTime=0.002;
filename='OriAnimation.gif';
h=figure;

for j=1:N
    clf   
    hold on
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    axis([-0.5*Vsize Vsize -0.5*Vsize Vsize -0.5*Vsize Vsize])
    
    Node=[pt1his(j,:);pt2his(j,:);pt3his(j,:)];    
    panel=[1 2 3];
    patch('Vertices',Node,'Faces',panel,'FaceColor','yellow');  
    
    Node=[PointCord(j,:);pt1his(j,:);pt2his(j,:);pt3his(j,:)];
    for i=1:4
        index=i;
        quiver3(Node(index,1), Node(index,2), Node(index,3)...
            ,rate*forceHis(j,index,1),rate*forceHis(j,index,2),rate*forceHis(j,index,3),'Color','red');
    end
    scatter3(PointCord(j,1),PointCord(j,2),PointCord(j,3),'filled');      
    pause(pauseTime);
    
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if j == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
    end 
end




