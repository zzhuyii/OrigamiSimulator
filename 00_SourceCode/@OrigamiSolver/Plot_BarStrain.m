%% Plot bar strain

function Plot_BarStrain(obj,loadController)

    View1=obj.viewAngle1;
    View2=obj.viewAngle2;
    Vsize=obj.displayRange;
    Vratio=obj.displayRangeRatio;

    % Plot Dots
    figure
    hold on
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    set(gcf, 'color', 'white');
    set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])
    
    A=size(Vsize);
    if A(1)==1    
        axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
    else
        axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
    end

    B=size(obj.barArea);
    BarNum=B(1);

%     barSx=loadController.barSxHis(end,:);
    barSx=loadController.barSxHis(end,1:96*8);
    
    maxSx = max(barSx);
    minSx = min(barSx);
    
    for i=1:96*8
        if barSx(i)>3/5*maxSx
            plot3([obj.newNode(obj.barConnect(i,1),1);
                obj.newNode(obj.barConnect(i,2),1)],...
                [obj.newNode(obj.barConnect(i,1),2);
                obj.newNode(obj.barConnect(i,2),2)],...
                [obj.newNode(obj.barConnect(i,1),3);
                obj.newNode(obj.barConnect(i,2),3)],...
                'Color','red','LineWidth',2);
        elseif barSx(i)>1/5*maxSx
            plot3([obj.newNode(obj.barConnect(i,1),1);
                obj.newNode(obj.barConnect(i,2),1)],...
                [obj.newNode(obj.barConnect(i,1),2);
                obj.newNode(obj.barConnect(i,2),2)],...
                [obj.newNode(obj.barConnect(i,1),3);
                obj.newNode(obj.barConnect(i,2),3)],...
                'Color',[1,0.7,0],'LineWidth',2);
        elseif barSx(i)>1/5*minSx
            plot3([obj.newNode(obj.barConnect(i,1),1);
                obj.newNode(obj.barConnect(i,2),1)],...
                [obj.newNode(obj.barConnect(i,1),2);
                obj.newNode(obj.barConnect(i,2),2)],...
                [obj.newNode(obj.barConnect(i,1),3);
                obj.newNode(obj.barConnect(i,2),3)],...
                'Color',[0.5,0.5,0.5],'LineWidth',2);
        elseif barSx(i)>3/5*minSx
            plot3([obj.newNode(obj.barConnect(i,1),1);
                obj.newNode(obj.barConnect(i,2),1)],...
                [obj.newNode(obj.barConnect(i,1),2);
                obj.newNode(obj.barConnect(i,2),2)],...
                [obj.newNode(obj.barConnect(i,1),3);
                obj.newNode(obj.barConnect(i,2),3)],...
                'Color',[0.5,0.7,1],'LineWidth',2);
        else
            plot3([obj.newNode(obj.barConnect(i,1),1);
                obj.newNode(obj.barConnect(i,2),1)],...
                [obj.newNode(obj.barConnect(i,1),2);
                obj.newNode(obj.barConnect(i,2),2)],...
                [obj.newNode(obj.barConnect(i,1),3);
                obj.newNode(obj.barConnect(i,2),3)],...
                'Color','blue','LineWidth',2);
        end
    end
    
    h=scatter3(-100,-100,-100,...
        'o','red','MarkerFaceColor','red');    
    
    h2=scatter3(-101,-100,-100,...
        'o','MarkerEdgeColor',[1,0.7,0], 'MarkerFaceColor',[1,0.7,0]);
    
    h3=scatter3(-102,-100,-100,...
        'o','yellow','MarkerFaceColor',[0.5,0.5,0.5]);    

    h4=scatter3(-103,-100,-100,...
        'o','cyan','MarkerFaceColor',[0.5,0.7,1]);
    
    h5=scatter3(-104,-100,-100,...
        'o','blue','MarkerFaceColor','blue');
    
    a=compose('%E  to %E',3/5*maxSx,maxSx);
    a2=compose('%E  to %E',1/5*maxSx,3/5*maxSx);
    a3=compose('%E  to %E',1/5*minSx,1/5*maxSx);
    a4=compose('%E  to %E',3/5*minSx,1/5*minSx);
    a5=compose('%E  to %E',minSx,3/5*minSx);
    
    legend([h,h2,h3,h4,h5],[a,a2,a3,a4,a5]);

    hold off
end
