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

function Plot_ContactForce(obj,undeformedNode,deformNode,contactForce)

    View1=obj.viewAngle1;
    View2=obj.viewAngle2;
    Vsize=obj.displayRange;
    Vratio=obj.displayRangeRatio;
    
    % Plot Dots
    figure
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    axis([-Vsize*Vratio Vsize -Vsize*Vratio Vsize -Vsize*Vratio Vsize])

    B=size(obj.newPanel);
    FaceNum=B(2);

    hold on
    for i=1:FaceNum
        tempPanel=cell2mat(obj.newPanel(i));
        patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor','yellow');
    end

    for i=1:FaceNum
        tempPanel=cell2mat(obj.newPanel(i));
        patch('Vertices',undeformedNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
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