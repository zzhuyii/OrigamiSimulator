%% Plot the loading and support reaction for the final step of loading

function plotLoadAndReaction(ViewControl,newNode,deformNode,newPanel,Load,Supp,NodeForce,LoadForce,PanelNum)

    View1=ViewControl(1);
    View2=ViewControl(2);
    Vsize=ViewControl(3);
    Vratio=ViewControl(4);
    
    % Plot Dots
    figure
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    axis([-Vsize*Vratio Vsize -Vsize*Vratio Vsize -Vsize*Vratio Vsize])

    B=size(newPanel);
    FaceNum=B(2);

    hold on
    for i=1:FaceNum
        tempPanel=cell2mat(newPanel(i));
        patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor','yellow');
    end

    for i=1:FaceNum
        tempPanel=cell2mat(newPanel(i));
        patch('Vertices',newNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
    end



    A=size(Load);
    loadNum=A(1);

    maxF=0;
    for i=1:loadNum
        index=Load(i,1);
        if norm(LoadForce(index,:))>maxF
            maxF=norm(LoadForce(index,:));
        end 
    end

    for i=1:loadNum
        index=Load(i,1);
        rate=norm(NodeForce(index,:))/maxF*Vsize;   
        quiver3(deformNode(index,1)-rate*NodeForce(index,1) ... 
        ,deformNode(index,2)-rate*NodeForce(index,2) ...
        ,deformNode(index,3)-rate*NodeForce(index,3)...
        ,rate*NodeForce(index,1),rate*NodeForce(index,2),rate*NodeForce(index,3),'Color','blue','AutoScale','off');
    end

    for i=1:loadNum
        index=Load(i,1);
        rate=10*norm(LoadForce(index,:))/maxF*Vsize;  
        quiver3(deformNode(index,1)-rate*LoadForce(index,1) ... 
        ,deformNode(index,2)-rate*LoadForce(index,2) ...
        ,deformNode(index,3)-rate*LoadForce(index,3)...
        ,rate*LoadForce(index,1),rate*LoadForce(index,2),rate*LoadForce(index,3),'Color','green','AutoScale','off');
    end

    A=size(Supp);
    supNum=A(1);

    for i=1:supNum
        index=Supp(i,1);
        rate=10*norm(NodeForce(index,:))/maxF*Vsize;  
        quiver3(deformNode(index,1)-rate*NodeForce(index,1) ... 
            ,deformNode(index,2)-rate*NodeForce(index,2) ...
            ,deformNode(index,3)-rate*NodeForce(index,3)...
            ,rate*NodeForce(index,1),rate*NodeForce(index,2),rate*NodeForce(index,3),'Color','red','AutoScale','off');
    end

    hold off
end