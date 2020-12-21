%% Calculate the bar length 
%
% Input:
%       newNode: the nodal coordinates of nodes after meshing;
%       barConnect: the two nodal number of each bars;
% Output: 
%       barLength: length of each bars;
%

function Mesh_BarLength(obj)
    A=size(obj.barConnect);
    barNum=A(1);
    obj.barLength=zeros(5,1);
    for i=1:barNum
        obj.barLength(i)=norm(obj.newNode(obj.barConnect(i,1),:)-...
            obj.newNode(obj.barConnect(i,2),:));
    end
end