%% plot the deformation history of the simulated results.

function plotDeformedHis(ViewControl,newNode,newPanel,UloadHis,UassembleHis,PanelNum)

View1=ViewControl(1);
View2=ViewControl(2);
Vsize=ViewControl(3);
Vratio=ViewControl(4);

pauseTime=0.05;
filename='OriAnimation.gif';

h=figure;
% view(View1,View2); 
% set(gca,'DataAspectRatio',[1 1 1])
% axis([-0.2*Vsize Vsize -0.2*Vsize Vsize -0.2*Vsize Vsize])
% axis([-Vsize Vsize -Vsize Vsize -Vsize Vsize])

B=size(newPanel);
FaceNum=B(2);

A=size(UassembleHis);
B=size(UloadHis);
Incre1=A(1);
Incre2=B(1);
n1=A(2);
n2=A(3);

set(gcf, 'color', 'white');
set(gcf,'Position',[0,-300,1000,1000])
    
for i=1:Incre1
    clf
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
    tempU=zeros(n1,n2);
    for j=1:n1
        for k=1:n2
           tempU(j,k)=UassembleHis(i,j,k);
        end
    end
    deformNode=newNode+tempU;
    for j=1:FaceNum
        tempPanel=cell2mat(newPanel(j));
        patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor','yellow');    
    end 
    pause(pauseTime); 

    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
    end 
end

for i=1:Incre2
    clf
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
    tempU=zeros(n1,n2);
    for j=1:n1
        for k=1:n2
           tempU(j,k)=UloadHis(i,j,k);
        end
    end
    deformNode=newNode+tempU;
    for j=1:FaceNum
        tempPanel=cell2mat(newPanel(j));
        patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor','yellow');     
    end
    pause(pauseTime);  
    
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
end
