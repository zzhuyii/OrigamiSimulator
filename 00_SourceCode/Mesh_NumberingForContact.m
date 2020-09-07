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

function [creaseRef]= Mesh_NumberingForContact(newCrease2OldCrease,oldCreaseNum)
    A=size(newCrease2OldCrease);
    Num=A(1);
    creaseRef=zeros(oldCreaseNum,8);
    for i=1:oldCreaseNum
        step=1;
        for j=1:Num
            if newCrease2OldCrease(j)==i 
                creaseRef(i,step)=j;
                step=step+1;
                if step>8
                    break
                end
            end
        end 
    end
end