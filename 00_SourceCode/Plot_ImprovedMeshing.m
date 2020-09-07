%% Plot the generated meshing of origami
%
% This funciton plots the origami after the meshing is generated. It is
% used for inspection and check the numbering system;
%
% Input:
%       viewControl: general control for plotting;
%       newNode: nodal coordinates of the generated mesh;
%       newPanel: panel setup for the new mesh;
%       barConnect: connectivity of bar elements
%

function Plot_ImprovedMeshing(viewControl,newNode,newPanel,barConnect)

    View1=viewControl(1);
    View2=viewControl(2);
    Vsize=viewControl(3);
    Vratio=viewControl(4);

    figure
    H=plot3(newNode(:,1)',newNode(:,2)',newNode(:,3)','.');
    hold on;
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize]);
    % Number Dots
    A=size(newNode);
    N=A(1);
    Sequence{1}='0';
    for i=1:N
        Sequence{i}=num2str(i);
    end
    text(newNode(:,1)',newNode(:,2)',newNode(:,3)',Sequence);
    A=size(barConnect);
    N=A(1);
    for i=1:N
        x=0.5*(newNode(barConnect(i,1),1)+newNode(barConnect(i,2),1));
        y=0.5*(newNode(barConnect(i,1),2)+newNode(barConnect(i,2),2));
        z=0.5*(newNode(barConnect(i,1),3)+newNode(barConnect(i,2),3));
        text(x,y,z,num2str(i),'Color','blue');
    end
    % Plot Panels
    B=size(newPanel);
    FaceNum=B(2);
    for i=1:FaceNum
        tempPanel=cell2mat(newPanel(i));
        patch('Vertices',newNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
    end
end