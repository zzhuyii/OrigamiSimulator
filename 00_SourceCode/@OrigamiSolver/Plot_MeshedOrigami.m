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

function Plot_MeshedOrigami(obj)

    figure
    H=plot3(obj.newNode(:,1)',obj.newNode(:,2)',obj.newNode(:,3)','.');
    hold on;
    view(obj.viewAngle1,obj.viewAngle2); 
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-obj.displayRangeRatio*obj.displayRange ... 
        obj.displayRange -obj.displayRangeRatio*obj.displayRange ... 
        obj.displayRange -obj.displayRangeRatio*obj.displayRange ...
        obj.displayRange]);
    % Number Dots
    A=size(obj.newNode);
    N=A(1);
    Sequence{1}='0';
    for i=1:N
        Sequence{i}=num2str(i);
    end
    if obj.showNumber==1
        text(obj.newNode(:,1)',obj.newNode(:,2)',obj.newNode(:,3)',Sequence);
    end
    A=size(obj.barConnect);
    N=A(1);
    if obj.showNumber==1
        for i=1:N
            x=0.5*(obj.newNode(obj.barConnect(i,1),1)+obj.newNode(obj.barConnect(i,2),1));
            y=0.5*(obj.newNode(obj.barConnect(i,1),2)+obj.newNode(obj.barConnect(i,2),2));
            z=0.5*(obj.newNode(obj.barConnect(i,1),3)+obj.newNode(obj.barConnect(i,2),3));
            text(x,y,z,num2str(i),'Color','blue');
        end
    end
    % Plot Panels
    B=size(obj.newPanel);
    FaceNum=B(2);
    for i=1:FaceNum
        tempPanel=cell2mat(obj.newPanel(i));
        patch('Vertices',obj.newNode,'Faces',tempPanel,...
            'EdgeColor',[0.5 0.5 0.5],...
            'FaceColor',obj.faceColorNumbering, ...
            'FaceAlpha',obj.faceAlphaNumbering);
    end
end