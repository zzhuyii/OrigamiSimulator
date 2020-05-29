%% Plot the locking force of the final steps

function plotLockForce(ViewControl,newNode,deformNode,newPanel,lockForce,PanelNum)

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

    A=size(lockForce);
    supNum=A(1);

    maxLock=0;
    for i=1:supNum
        if norm(lockForce(i,:))>maxLock
            maxLock=norm(lockForce(i,:));
        end    
    end
    
    for i=1:supNum
        index=i;    
        rate=norm(lockForce(index,:))/maxLock*10*Vsize;    
        quiver3(deformNode(index,1) ... 
            ,deformNode(index,2) ...
            ,deformNode(index,3)...
            ,rate*lockForce(index,1),rate*lockForce(index,2),rate*lockForce(index,3),'Color','red','AutoScale','off');
    end
    
    hold off
end