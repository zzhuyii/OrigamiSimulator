%% This function adds a bar to the origami model

function AddHinge(obj,ni,nj,nk,nl,BarArea,L,SprStiff,RotZero)

    currentBarNum=size(obj.barType);
    currentBarNum=max(currentBarNum);

    obj.barType(currentBarNum+1,:)=1;
    obj.barConnect(currentBarNum+1,:)=[nj,nk];
    obj.barArea(currentBarNum+1,:)=BarArea;
    obj.barLength(currentBarNum+1,:)=L;

    obj.sprK(currentBarNum+1,:)=SprStiff;
    obj.sprIJKL(currentBarNum+1,:)=[ni,nj,nk,nl];
    obj.currentSprZeroStrain(currentBarNum+1,:)=RotZero;

end