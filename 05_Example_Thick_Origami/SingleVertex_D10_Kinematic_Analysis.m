%% Single Vertex Kinematical Animation
clear all
close all
clc

% This is the value we used to generate the Yoshimura vertex
fold1=0.2;
fold2=0.245;

fold1=0.9;
fold2=0.92;

fold1=0;
fold2=0;

sectorAngle=[0.125*pi,0.5*pi,0.125*pi,0.125*pi,0.125*pi,0.125*pi,0.5*pi,0.125*pi,0.125*pi,0.125*pi];
foldAngle=[-fold1*pi,fold2*pi,fold2*pi,-fold1*pi,fold1*pi,-fold1*pi,fold2*pi,fold2*pi,-fold1*pi,fold1*pi];
t=[0.2,0,-0.2,0.2,-0.2,0.2,0,-0.2,0.2,-0.2];

draw(sectorAngle,foldAngle,t)




%% This function computee the derivative of constraint function
% Because of the structure of D-H transformation matrix, we just need to
% compute the top 3 rows
function D=DiffConstraintFunc(sectorAngle,foldAngle,a,h)

    T=ConstraintFunc(sectorAngle,foldAngle,a);
    
    D=zeros(12,10); 
    % Numerical derivation of the constraint matrix this is a 12 by 10 
    % matrix, as we have 12 constraint equation and 10 variable

    for j=1:10
        tempFoldAngle1=foldAngle;
        tempFoldAngle2=foldAngle;

        tempFoldAngle1(j)=foldAngle(j)+h;
        tempFoldAngle2(j)=foldAngle(j)-h;

        T1=ConstraintFunc(sectorAngle,tempFoldAngle1,a);
        T2=ConstraintFunc(sectorAngle,tempFoldAngle2,a);

        Dmat=(T2(1:3,:)-T1(1:3,:))/2/h;

        Dmat=reshape(Dmat,12,1);
        D(:,j)=Dmat;

    end

end

%% This function computes the constraint function's value
% When the output value gives an zero matrix, we know that the constriant
% is met

function T=ConstraintFunc(sectorAngle,foldAngle,a)
    tempT=[1 0 0 0;
           0 1 0 0;
           0 0 1 0;
           0 0 0 1];
    
    for i=1:10
        tempT=tempT*Tmat(sectorAngle(i),foldAngle(i),a(i),0); 
        % This is the main crease vector 
    end
    
    T=tempT-eye(4);
end



%% D-H transformation function
function T=Tmat(alpha,theta,a,d)

    T=[cos(theta)  -sin(theta)*cos(alpha)  sin(theta)*sin(alpha)    a*cos(theta);
       sin(theta)  cos(theta)*cos(alpha)   -cos(theta)*sin(alpha)   a*sin(theta);
       0           sin(alpha)              cos(alpha)               d;
       0           0                       0                        1];
    
end


%% Plot the geometry of vertex
function draw(sectorAngle,foldAngle,a)

    e0=zeros(3,3);
    e1=zeros(3,3);
    e2=zeros(3,3);

    e3=zeros(3,3);
    e4=zeros(3,3);
    e5=zeros(3,3);
    
    tempT=[1 0 0 0;
           0 1 0 0;
           0 0 1 0;
           0 0 0 1];
    
    tempT1=tempT*Tmat(0,0,0,1);
    tempT2=tempT*Tmat(sectorAngle(1),0,0,0)*Tmat(0,0,0,1);
    
    e1(:,1)=tempT1(1:3,4);
    e2(:,1)=tempT2(1:3,4);
    
    figure
    hold on
    view(135,45)
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-1.5, 1.5, -1.5, 1.5, -1.5, 1.5]);

    for i=1:10
    
        avalue=a(1);
        if i==2
            tempT3=tempT*Tmat(0,foldAngle(i),-avalue,0)*Tmat(0,0,0,1);
            tempT4=tempT*Tmat(sectorAngle(i),foldAngle(i),-avalue,0)*Tmat(0,0,0,1);
            tempT5=tempT*Tmat(0,foldAngle(i),-avalue,0);
        elseif i==7
            tempT3=tempT*Tmat(0,foldAngle(i),-avalue,0)*Tmat(0,0,0,1);
            tempT4=tempT*Tmat(sectorAngle(i),foldAngle(i),-avalue,0)*Tmat(0,0,0,1);
            tempT5=tempT*Tmat(0,foldAngle(i),-avalue,0);
        else
            tempT3=tempT*Tmat(0,foldAngle(i),0,0)*Tmat(0,0,0,1);
            tempT4=tempT*Tmat(sectorAngle(i),foldAngle(i),0,0)*Tmat(0,0,0,1);
            tempT5=tempT;
        end
    

        tempT2=tempT*Tmat(0,foldAngle(i),a(i),0)*Tmat(0,0,0,1);
        tempT=tempT*Tmat(sectorAngle(i),foldAngle(i),a(i),0);      
        tempT1=tempT*Tmat(0,0,0,1);    

    
        % node of one facet panel
        e0(:,i+1)=tempT(1:3,4);
        e1(:,i+1)=tempT1(1:3,4);
        e2(:,i+1)=tempT2(1:3,4);
    
        % node of another facet panel
        e3(:,i+1)=tempT3(1:3,4);
        e4(:,i+1)=tempT4(1:3,4);
        e5(:,i+1)=tempT5(1:3,4);
    
    
        % Facet Panels
        patch([e0(1,i+1),e1(1,i+1),e2(1,i+1)], ...
          [e0(2,i+1),e1(2,i+1),e2(2,i+1)], ...
          [e0(3,i+1),e1(3,i+1),e2(3,i+1)],'yellow')
    
        patch([e3(1,i+1),e4(1,i+1),e5(1,i+1)], ...
          [e3(2,i+1),e4(2,i+1),e5(2,i+1)], ...
          [e3(3,i+1),e4(3,i+1),e5(3,i+1)],'yellow')


        % Side Panels
        patch([e0(1,i+1),e2(1,i+1),e3(1,i+1),e5(1,i+1)], ...
          [e0(2,i+1),e2(2,i+1),e3(2,i+1),e5(2,i+1)], ...
          [e0(3,i+1),e2(3,i+1),e3(3,i+1),e5(3,i+1)],'yellow')

        % Side Panels
        patch([e1(1,i+1),e2(1,i+1),e3(1,i+1),e4(1,i+1)], ...
          [e1(2,i+1),e2(2,i+1),e3(2,i+1),e4(2,i+1)], ...
          [e1(3,i+1),e2(3,i+1),e3(3,i+1),e4(3,i+1)],'yellow')

        % Side Panels
        patch([e0(1,i+1),e1(1,i+1),e4(1,i+1),e5(1,i+1)], ...
          [e0(2,i+1),e1(2,i+1),e4(2,i+1),e5(2,i+1)], ...
          [e0(3,i+1),e1(3,i+1),e4(3,i+1),e5(3,i+1)],'yellow')
    
    end

end