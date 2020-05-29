%% Plot the original pattern design before the automated meshing

function plotOriginalMeshing(Node,Panel,CreaseNum,Crease,ViewControl)
    View1=ViewControl(1);
    View2=ViewControl(2);
    Vsize=ViewControl(3);
    Vratio=ViewControl(4);

    figure
    H=plot3(Node(:,1)',Node(:,2)',Node(:,3)','.');
    hold on;
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-Vsize*Vratio Vsize -Vsize*Vratio Vsize -Vsize*Vratio Vsize]);
    % Number Dots
    A=size(Node);
    N=A(1);
    Sequence{1}='0';
    for i=1:N
        Sequence{i}=num2str(i);
    end
    text(Node(:,1)',Node(:,2)',Node(:,3)',Sequence);
    for i=1:CreaseNum
        x=0.5*(Node(Crease(i,1),1)+Node(Crease(i,2),1));
        y=0.5*(Node(Crease(i,1),2)+Node(Crease(i,2),2));
        z=0.5*(Node(Crease(i,1),3)+Node(Crease(i,2),3));
        text(x,y,z,num2str(i),'Color','blue');
    end
    % Plot Panels
    B=size(Panel);
    FaceNum=B(2);
    for i=1:FaceNum
        tempPanel=cell2mat(Panel(i));
        patch('Vertices',Node,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
    end
end