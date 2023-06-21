%% Plot the deformation history and temperature changes
%
% This code plot the deformation history of the origami system and overlays
% the changes of nodal temperatures on top.
%
% Input:
%       viewControl: general control input of plotting;
%       newNode: the nodal coordinates of origami;
%       UhisThermal: deformation history;
%       newPanel: panel connectivity;
%       TemperatureHistory: temperature history of nodes;
%

function Plot_DeformedHisTemp(obj,beforeLoadingNode,UhisThermal,temperatureHistory)

    View1=obj.viewAngle1;
    View2=obj.viewAngle2;
    Vsize=obj.displayRange;
    Vratio=obj.displayRangeRatio;
        
    A=size(UhisThermal);
    ThermalStep=A(1);
    nodeNum=A(2);
    
    maxT = max(temperatureHistory(:));
    minT = min(temperatureHistory(:));
    
    pauseTime=obj.holdTime;
    filename='OriThermalAnimation.gif';
    
    MyFig=figure;
    B=size(obj.newPanel);
    FaceNum=B(2);    

    set(gcf, 'color', 'white');
    set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])
        
    for step=1:ThermalStep
        
        clf
        hold on
        view(View1,View2); 

        A=size(Vsize);
        if A(1)==1    
            axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
        else
            axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
        end
    
        set(gca,'DataAspectRatio',[1 1 1])
        deformNode=beforeLoadingNode+squeeze(UhisThermal(step,:,:));
        for i=1:FaceNum
            tempPanel=cell2mat(obj.newPanel(i));
            patch('Vertices',deformNode,'Faces',tempPanel, ...
            'FaceColor',obj.faceColorAnimation, ...
            'FaceAlpha', obj.faceAlphaAnimation);
        end

        for i=1:FaceNum
            tempPanel=cell2mat(obj.newPanel(i));
            patch('Vertices',beforeLoadingNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
        end
        
        T=squeeze(temperatureHistory(:,step));
        
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

        h1=scatter3(-100,-100,-100,...
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
        
        legend([h1,h2,h3,h4,h5],[a,a2,a3,a4,a5]);
        
        hold off
        pause(pauseTime); 

        frame = getframe(MyFig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File 
        if step == 1 
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime', pauseTime); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
        end 
    end   
end
