%% Plot the deformed configuration with temperature of node
%
% This code plot the deformed configuration and overlay the temperature of
% the node onto the deformed configuration
%
% Input:
%       viewControl: general control input of plotting;
%       newNode: the nodal coordinates of origami;
%       deformNode: the deformed nodal coordinates;
%       newPanel: panel connectivity;
%       T: temperature profile of nodes;
%       thermalBoundaryPanel: the panels set as thermal boundary;
%       newPanel2OldPanel: mapping between new and old panel numbers;
%

function Plot_DeformedShapeTemp(obj,thermal,newNode,deformNode,T)

    View1=obj.viewAngle1;
    View2=obj.viewAngle2;
    Vsize=obj.displayRange;
    Vratio=obj.displayRangeRatio;
    
    A=size(thermal.thermalBoundaryPanelVec);
    BCnum=A(1);
        
    % Plot Dots
    figure
    hold on
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    axis([-Vsize*Vratio Vsize -Vsize*Vratio Vsize -Vsize*Vratio Vsize])

    B=size(obj.newPanel);
    FaceNum=B(2);

    for i=1:FaceNum
        tempPanel=cell2mat(obj.newPanel(i));
        

        check=0;
        for j=1:BCnum
            if obj.newPanel2OldPanel(i)==...
                    thermal.thermalBoundaryPanelVec(j) && ...
                    obj.newPanel2OldPanel(i)~=0
                check=1;
            end
        end
        if check==1
            patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor',[0.5 0.7 0.7]);
        else
            patch('Vertices',deformNode,'Faces',tempPanel,...
                'FaceColor',obj.faceColorAnimation,...
                'FaceAlpha', obj.faceAlphaAnimation);
        end
    end

    for i=1:FaceNum
        tempPanel=cell2mat(obj.newPanel(i));
        patch('Vertices',newNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
    end

    A=size(T);
    nodeNum=A(1);
    
    maxT = max(T);
    minT = min(T);
    
    for i=1:nodeNum
        if T(i)>4/5*(maxT-minT)+minT
            scatter3(deformNode(i,1),deformNode(i,2),deformNode(i,3),...
                'o','red','MarkerFaceColor','red');
        elseif T(i)>3/5*(maxT-minT)+minT
            scatter3(deformNode(i,1),deformNode(i,2),deformNode(i,3),...
                'o','MarkerEdgeColor',[1,0.7,0], 'MarkerFaceColor',[1,0.7,0])    
        elseif T(i)>2/5*(maxT-minT)+minT
            scatter3(deformNode(i,1),deformNode(i,2),deformNode(i,3),...
                'o','yellow','MarkerFaceColor','yellow') 
        elseif T(i)>1/5*(maxT-minT)+minT
            scatter3(deformNode(i,1),deformNode(i,2),deformNode(i,3),...
                'o','cyan','MarkerFaceColor','cyan') 
        else
            scatter3(deformNode(i,1),deformNode(i,2),deformNode(i,3),...
                'o','blue','MarkerFaceColor','blue') 
        end
    end
    
    h=scatter3(-100,-100,-100,...
        'o','red','MarkerFaceColor','red');    
    
    h2=scatter3(-101,-100,-100,...
        'o','MarkerEdgeColor',[1,0.7,0], 'MarkerFaceColor',[1,0.7,0]);
    
    h3=scatter3(-102,-100,-100,...
        'o','yellow','MarkerFaceColor','yellow');    

    h4=scatter3(-103,-100,-100,...
        'o','cyan','MarkerFaceColor','cyan');
    
    h5=scatter3(-104,-100,-100,...
        'o','blue','MarkerFaceColor','blue');
    
    a=compose('%9.2f  to %9.2f',4/5*(maxT-minT)+minT,maxT);
    a2=compose('%9.2f  to %9.2f',3/5*(maxT-minT)+minT,4/5*(maxT-minT)+minT);
    a3=compose('%9.2f  to %9.2f',2/5*(maxT-minT)+minT,3/5*(maxT-minT)+minT);
    a4=compose('%9.2f  to %9.2f',1/5*(maxT-minT)+minT,2/5*(maxT-minT)+minT);
    a5=compose('%9.2f  to %9.2f',0/5*(maxT-minT)+minT,1/5*(maxT-minT)+minT);
    
    legend([h,h2,h3,h4,h5],[a,a2,a3,a4,a5]);

    hold off
end
