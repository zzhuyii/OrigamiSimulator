%% Plot the locking force of the final steps
%
% The function plots the deformed origami shape and overlay the contact
% forces on top of the figure with arrows. 
%
% Input:
%       viewControl: general plotting set up
%       newNode: nodal coordinates with unloaded origami
%       deformNode: deformed nodal coordinates;
%       newPanel: panel information;
%       contactForce: force vector for contact;
%

function Plot_ContactForce(viewControl,newNode,deformNode,...
    newPanel,contactForce)

    View1=viewControl(1);
    View2=viewControl(2);
    Vsize=viewControl(3);
    Vratio=viewControl(4);
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

    A=size(contactForce);
    supNum=A(1);

    maxLock=0;
    for i=1:supNum
        if norm(contactForce(i,:))>maxLock
            maxLock=norm(contactForce(i,:));
        end    
    end
    
    for i=1:supNum
        index=i;    
        rate=norm(contactForce(index,:))/maxLock*10*Vsize;    
        quiver3(deformNode(index,1) ... 
            ,deformNode(index,2) ...
            ,deformNode(index,3)...
            ,rate*contactForce(index,1),rate*contactForce(index,2),rate*contactForce(index,3),'Color','red','AutoScale','off');
    end
    
    hold off
end