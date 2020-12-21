%% Plot the loading and support reaction for the final step of loading
%
% This function plots the support reaction and the loading forces acting
% on the origami structure.
%
% Input:
%       viewControl: general plotting set up;
%       newNode: the unloaded nodal coordinates;
%       deformNode: the deformed nodal coordinates;
%       newPanel: panel information;
%       load: load vector when setting up loading;
%       supp: support information;
%       nodeForce: nodal force vector;
%       loadForce: load force vector;
%

function Plot_LoadAndReaction(obj,undeformedNode,deformNode,...
    load,supp,nodeForce,loadForce)

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



    A=size(load);
    loadNum=A(1);

    maxF=0;
    for i=1:loadNum
        index=load(i,1);
        if norm(loadForce(index,:))>maxF
            maxF=norm(loadForce(index,:));
        end 
    end

    for i=1:loadNum
        index=load(i,1);
        rate=norm(nodeForce(index,:))/maxF*Vsize;   
        quiver3(deformNode(index,1)-rate*nodeForce(index,1) ... 
        ,deformNode(index,2)-rate*nodeForce(index,2) ...
        ,deformNode(index,3)-rate*nodeForce(index,3)...
        ,rate*nodeForce(index,1),rate*nodeForce(index,2),rate*nodeForce(index,3),'Color','blue','AutoScale','off');
    end

    for i=1:loadNum
        index=load(i,1);
        rate=10*norm(loadForce(index,:))/maxF*Vsize;  
        quiver3(deformNode(index,1)-rate*loadForce(index,1) ... 
        ,deformNode(index,2)-rate*loadForce(index,2) ...
        ,deformNode(index,3)-rate*loadForce(index,3)...
        ,rate*loadForce(index,1),rate*loadForce(index,2),rate*loadForce(index,3),'Color','green','AutoScale','off');
    end

    A=size(supp);
    supNum=A(1);

    for i=1:supNum
        index=supp(i,1);
        rate=10*norm(nodeForce(index,:))/maxF*Vsize;  
        quiver3(deformNode(index,1)-rate*nodeForce(index,1) ... 
            ,deformNode(index,2)-rate*nodeForce(index,2) ...
            ,deformNode(index,3)-rate*nodeForce(index,3)...
            ,rate*nodeForce(index,1),rate*nodeForce(index,2),rate*nodeForce(index,3),'Color','red','AutoScale','off');
    end

    hold off
end