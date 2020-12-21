%% Generates a number reference info for calculating the locking
%
% Input include:
%    {newCrease2OldCrease}: how the new creases(bars) are connectd; 
%    [oldCreaseNum]: total number of old creases;
%
% Output include: 
%    [creaseRef]: how ith old crease are further divided into different
%                 new rotaional springs;
%

function Mesh_NumberingForContact(obj)
    A=size(obj.newCrease2OldCrease);
    Num=A(1);
    obj.creaseRef=zeros(obj.oldCreaseNum,8);
    for i=1:obj.oldCreaseNum
        step=1;
        for j=1:Num
            if obj.newCrease2OldCrease(j)==i 
                obj.creaseRef(i,step)=j;
                step=step+1;
                if step>8
                    break
                end
            end
        end 
    end
end