%% Calculate the bar length 
%
% Input:
%       newNode: the nodal coordinates of nodes after meshing;
%       barConnect: the two nodal number of each bars;
% Output: 
%       barLength: length of each bars;
%

function barLength=Mesh_BarLength(newNode,barConnect)
    A=size(barConnect);
    barNum=A(1);
    barLength=zeros(5,1);
    for i=1:barNum
        barLength(i)=norm(newNode(barConnect(i,1),:)-...
            newNode(barConnect(i,2),:));
    end
end