%% Plot the deformed configuration
%
% This function plots the deformed shape of the origami without plotting
% the undeformed (refernece) origami configuration.
%
% Input:
%       viewControl: general view control setup;
%       deformNode: the deformed nodal coordinates;
%       newPanel: panel information;
%

function Plot_DeformedShapeOnly(viewControl,deformNode,newPanel)

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
