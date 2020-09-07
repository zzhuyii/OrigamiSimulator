%% Plot the deformed configuration of the simulation
%
% This code plots the deformed confirguration of an origami system and the
% reference configuration of the origami before loading). 
%
% Input:
%       viewControl: general plotting setup;
%       newNode: reference nodal coordinates;
%       deformNode: deformed nodal coordinates;
%       newPanel: panel information;
%

function Plot_DeformedShape(viewControl,newNode,deformNode,newPanel)

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

for i=1:FaceNum
    tempPanel=cell2mat(newPanel(i));
    patch('Vertices',deformNode,'Faces',tempPanel,'FaceColor','yellow');
end

for i=1:FaceNum
    tempPanel=cell2mat(newPanel(i));
    patch('Vertices',newNode,'Faces',tempPanel,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
end