%% This function adds a bar to the origami model

function AddBar(obj,n1,n2,BarArea,L)

    currentBarNum=size(obj.barType);
    currentBarNum=max(currentBarNum);

    obj.barType(currentBarNum+1,:)=1;
    obj.barConnect(currentBarNum+1,:)=[n1,n2];
    obj.barArea(currentBarNum+1,:)=BarArea;
    obj.barLength(currentBarNum+1,:)=L;

    obj.sprK(currentBarNum+1,:)=0;
    obj.sprIJKL(currentBarNum+1,:)=[0,0,0,0];
    obj.currentSprZeroStrain(currentBarNum+1,:)=0;


end