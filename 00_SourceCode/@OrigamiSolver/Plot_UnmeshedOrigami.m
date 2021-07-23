%% Plot the original pattern design before generating the new meshing
%
% This funciton generates the figure of the input origami pattern before
% generating the new meshing. The code is used for inspection.
%
% Input:
%       node0: input nodal coordinates;
%       panel0: input panel setup;
%       oldCreaseNum: total number of creases;
%       oldCreaseConnect: crease connectivity;
%       viewControl: general plotting setup;
%

function Plot_UnmeshedOrigami(obj)
    figure
    H=plot3(obj.node0(:,1)',obj.node0(:,2)',obj.node0(:,3)','.');
    hold on;
    view(obj.viewAngle1,obj.viewAngle2); 
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-obj.displayRange*obj.displayRangeRatio ...
        obj.displayRange -obj.displayRange*obj.displayRangeRatio ...
        obj.displayRange -obj.displayRange*obj.displayRangeRatio ...
        obj.displayRange]);
    % Number Dots
    A=size(obj.node0);
    N=A(1);
    Sequence{1}='0';
    for i=1:N
        Sequence{i}=num2str(i);
    end
    
    if obj.showNumber==1
        text(obj.node0(:,1)',obj.node0(:,2)',obj.node0(:,3)',Sequence);
        for i=1:obj.oldCreaseNum
            x=0.5*(obj.node0(obj.oldCreaseConnect(i,1),1)+...
                obj.node0(obj.oldCreaseConnect(i,2),1));
            y=0.5*(obj.node0(obj.oldCreaseConnect(i,1),2)+...
                obj.node0(obj.oldCreaseConnect(i,2),2));
            z=0.5*(obj.node0(obj.oldCreaseConnect(i,1),3)+...
                obj.node0(obj.oldCreaseConnect(i,2),3));
            text(x,y,z,num2str(i),'Color','blue');
        end
    end
    % Plot Panels
    B=size(obj.panel0);
    FaceNum=B(2);
    for i=1:FaceNum
        tempPanel=cell2mat(obj.panel0(i));
        patch('Vertices',obj.node0,'Faces',...
            tempPanel,'EdgeColor',[0.5 0.5 0.5],...
            'FaceColor',obj.faceColorNumbering, ...
            'FaceAlpha',obj.faceAlphaNumbering);
    end
end