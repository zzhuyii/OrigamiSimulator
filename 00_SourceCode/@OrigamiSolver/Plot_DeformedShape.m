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

function Plot_DeformedShape(obj,undeformedNode,deformNode)

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

for i=1:FaceNum
    tempPanel=cell2mat(obj.newPanel(i));
    patch('Vertices',deformNode,'Faces',tempPanel,...
        'FaceColor',obj.faceColorAnimation, ...
        'FaceAlpha',obj.faceAlphaAnimation);
end

for i=1:FaceNum
    tempPanel=cell2mat(obj.newPanel(i));
    patch('Vertices',undeformedNode,'Faces',tempPanel,...
        'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0);
end