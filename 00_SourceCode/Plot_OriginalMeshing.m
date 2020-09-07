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

function Plot_OriginalMeshing(node0,panel0,oldCreaseNum,oldCreaseConnect,viewControl)
    View1=viewControl(1);
    View2=viewControl(2);
    Vsize=viewControl(3);
    Vratio=viewControl(4);

    figure
    H=plot3(node0(:,1)',node0(:,2)',node0(:,3)','.');
    hold on;
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-Vsize*Vratio Vsize -Vsize*Vratio Vsize -Vsize*Vratio Vsize]);
    % Number Dots
    A=size(node0);
    N=A(1);
    Sequence{1}='0';
    for i=1:N
        Sequence{i}=num2str(i);
    end
    text(node0(:,1)',node0(:,2)',node0(:,3)',Sequence);
    for i=1:oldCreaseNum
        x=0.5*(node0(oldCreaseConnect(i,1),1)+node0(oldCreaseConnect(i,2),1));
        y=0.5*(node0(oldCreaseConnect(i,1),2)+node0(oldCreaseConnect(i,2),2));
        z=0.5*(node0(oldCreaseConnect(i,1),3)+node0(oldCreaseConnect(i,2),3));
        text(x,y,z,num2str(i),'Color','blue');
    end
    % Plot Panels
    B=size(panel0);
    FaceNum=B(2);
    for i=1:FaceNum
        tempPanel=cell2mat(panel0(i));
        patch('Vertices',node0,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
    end
end