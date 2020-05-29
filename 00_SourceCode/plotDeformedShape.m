%% Plot the deformed configuration of the simulation

function plotDeformedShape(ViewControl,newNode,deformNode,newPanel,PanelNum)

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

for i=1:FaceNum
    tempPanel=cell2mat(newPanel(i));
    patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor','yellow');
end

for i=1:FaceNum
    tempPanel=cell2mat(newPanel(i));
    patch('Vertices',newNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
end