%% Three-Point Spring Assemble Forces

function [TPsprForce,TPsprStiff] = ThreePointSpring_GlobalForceAndStiff(obj,currentU,newNode,sprK,sprNode,sprTheta0)

    A=size(sprK);
    TPsprNodeNum=A(1);

    A=size(currentU);
    NodeNum=A(1);
    
    TPsprForce=zeros(3*NodeNum,1);  
    TPsprStiff=zeros(3*NodeNum,3*NodeNum);

    for i=1:TPsprNodeNum
        node1Num=sprNode(i,1);
        node2Num=sprNode(i,2);
        node3Num=sprNode(i,3);

        tempF=findForce(sprK(i),sprTheta0(i),...
            newNode(node1Num,:)+currentU(node1Num,:), ...
            newNode(node2Num,:)+currentU(node2Num,:), ...
            newNode(node3Num,:)+currentU(node3Num,:));

        TPsprForce(3*node1Num-2:3*node1Num)=tempF(1:3);
        TPsprForce(3*node2Num-2:3*node2Num)=tempF(4:6);
        TPsprForce(3*node3Num-2:3*node3Num)=tempF(7:9);


        tempK=findStiff(sprK(i),sprTheta0(i),...
            newNode(node1Num,:)+currentU(node1Num,:), ...
            newNode(node2Num,:)+currentU(node2Num,:), ...
            newNode(node3Num,:)+currentU(node3Num,:));

        tempK=(tempK+tempK')/2;

        TPsprStiff(3*node1Num-2:3*node1Num,3*node1Num-2:3*node1Num)=tempK(1:3,1:3);
        TPsprStiff(3*node1Num-2:3*node1Num,3*node2Num-2:3*node2Num)=tempK(1:3,4:6);
        TPsprStiff(3*node1Num-2:3*node1Num,3*node3Num-2:3*node3Num)=tempK(1:3,7:9);
        TPsprStiff(3*node2Num-2:3*node2Num,3*node1Num-2:3*node1Num)=tempK(4:6,1:3);
        TPsprStiff(3*node2Num-2:3*node2Num,3*node2Num-2:3*node2Num)=tempK(4:6,4:6);
        TPsprStiff(3*node2Num-2:3*node2Num,3*node3Num-2:3*node3Num)=tempK(4:6,7:9);
        TPsprStiff(3*node3Num-2:3*node3Num,3*node1Num-2:3*node1Num)=tempK(7:9,1:3);
        TPsprStiff(3*node3Num-2:3*node3Num,3*node2Num-2:3*node2Num)=tempK(7:9,4:6);
        TPsprStiff(3*node3Num-2:3*node3Num,3*node3Num-2:3*node3Num)=tempK(7:9,7:9);
    end
end



function Uspr=SprPotential(sprK,sprTheta0,node11,node12,node13,...
    node21,node22,node23,node31,node32,node33)

    node1=[node11,node12,node13];
    node2=[node21,node22,node23];
    node3=[node31,node32,node33];

    v1=node1-node2;
    v2=node3-node2;

    v1=v1/norm(v1);
    v2=v2/norm(v2);

    theta=acos(dot(v1,v2));
    Uspr=0.5*sprK*(theta-sprTheta0)^2;
end


function Fspr=findForce(sprK,sprTheta0,node1,node2,node3)
    Fspr=zeros(9,1);
    delta=10^-8;


    Fspr(1)=(SprPotential(sprK,sprTheta0,node1(1)+delta,node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1)-delta,node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3)))/2/delta;

    Fspr(2)=(SprPotential(sprK,sprTheta0,node1(1),node1(2)+delta,node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2)-delta,node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3)))/2/delta;

    Fspr(3)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3)+delta,...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3)-delta,...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3)))/2/delta;

    Fspr(4)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1)+delta,node2(2),node2(3),node3(1),node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1)-delta,node2(2),node2(3),node3(1),node3(2),node3(3)))/delta;

    Fspr(5)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2)+delta,node2(3),node3(1),node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2)-delta,node2(3),node3(1),node3(2),node3(3)))/delta;

    Fspr(6)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3)+delta,node3(1),node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3)-delta,node3(1),node3(2),node3(3)))/delta;

    Fspr(7)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1)+delta,node3(2),node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1)-delta,node3(2),node3(3)))/delta;

    Fspr(8)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2)+delta,node3(3))-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2)-delta,node3(3)))/delta;

    Fspr(9)=(SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3)+delta)-...
        SprPotential(sprK,sprTheta0,node1(1),node1(2),node1(3),...
        node2(1),node2(2),node2(3),node3(1),node3(2),node3(3)-delta))/delta;
end


function stiffMat=findStiff(sprK,sprTheta0,node1,node2,node3)
    stiffMat=zeros(9);
    delta=10^-8;
    baseF=findForce(sprK,sprTheta0,node1,node2,node3);


    node1d1p=node1;
    node1d1n=node1;
    node1d1p(1)=node1(1)+delta;
    node1d1n(1)=node1(1)-delta;
    stiffMat(1,:)=(findForce(sprK,sprTheta0,node1d1p,node2,node3)-...
        findForce(sprK,sprTheta0,node1d1n,node2,node3))/2/delta;

    node1d2p=node1;
    node1d2p(2)=node1(2)+delta;
    node1d2n=node1;
    node1d2n(2)=node1(2)-delta;
    stiffMat(2,:)=(findForce(sprK,sprTheta0,node1d2p,node2,node3)-...
        findForce(sprK,sprTheta0,node1d2n,node2,node3))/2/delta;
    
    node1d3p=node1;
    node1d3p(3)=node1(3)+delta;
    node1d3n=node1;
    node1d3n(3)=node1(3)-delta;
    stiffMat(3,:)=(findForce(sprK,sprTheta0,node1d3p,node2,node3)-...
        findForce(sprK,sprTheta0,node1d3n,node2,node3))/2/delta;

    node2d1p=node2;
    node2d1p(1)=node2(1)+delta;
    node2d1n=node2;
    node2d1n(1)=node2(1)-delta;
    stiffMat(4,:)=(findForce(sprK,sprTheta0,node1,node2d1p,node3)-...
        findForce(sprK,sprTheta0,node1,node2d1n,node3))/2/delta;

    node2d2p=node2;
    node2d2p(2)=node2(2)+delta;
    node2d2n=node2;
    node2d2n(2)=node2(2)-delta;
    stiffMat(5,:)=(findForce(sprK,sprTheta0,node1,node2d2p,node3)-...
        findForce(sprK,sprTheta0,node1,node2d2n,node3))/2/delta;

    node2d3p=node2;
    node2d3p(3)=node2(3)+delta;
    node2d3n=node2;
    node2d3n(3)=node2(3)-delta;
    stiffMat(6,:)=(findForce(sprK,sprTheta0,node1,node2d3p,node3)-...
        findForce(sprK,sprTheta0,node1,node2d3n,node3))/2/delta;

    node3d1p=node3;
    node3d1p(1)=node3(1)+delta;
    node3d1n=node3;
    node3d1n(1)=node3(1)-delta;
    stiffMat(7,:)=(findForce(sprK,sprTheta0,node1,node2,node3d1p)-...
        findForce(sprK,sprTheta0,node1,node2,node3d1n))/2/delta;

    node3d2p=node3;
    node3d2p(2)=node3(2)+delta;
    node3d2n=node3;
    node3d2n(2)=node3(2)-delta;
    stiffMat(8,:)=(findForce(sprK,sprTheta0,node1,node2,node3d2p)-...
        findForce(sprK,sprTheta0,node1,node2,node3d2n))/2/delta;

    node3d3p=node3;
    node3d3p(3)=node3(3)+delta;
    node3d3n=node3;
    node3d3n(3)=node3(3)-delta;
    stiffMat(9,:)=(findForce(sprK,sprTheta0,node1,node2,node3d3p)-...
        findForce(sprK,sprTheta0,node1,node2,node3d3n))/2/delta;

%     stiffMat=(stiffMat+stiffMat')/2;
end