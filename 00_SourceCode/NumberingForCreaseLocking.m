%% Generates a number reference info for calculating the locking

function [CreaseRef]= NumberingForCreaseLocking(oldCrease,CreaseNum,BarType)
    A=size(oldCrease);
    Num=A(1);
    CreaseRef=zeros(CreaseNum,8);
    for i=1:CreaseNum
        step=1;
        for j=1:Num
            if oldCrease(j)==i 
                CreaseRef(i,step)=j;
                step=step+1;
                if step>8
                    break
                end
            end
        end 
    end
end