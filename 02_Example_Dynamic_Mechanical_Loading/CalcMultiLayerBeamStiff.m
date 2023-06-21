%% This code solve the EIcomp and the centroid of multi-layer creases

function[EIcomp,C]=CalcMultiLayerBeamStiff(thickness,E,fiberNum)

    totalThick=sum(thickness);    
    tempLayer=1;
    tempThick=thickness(1);    
    tempE=E(1);
    
    % First we solve for the centroid
    S=0;
    A=0;
    
    for i=1:fiberNum
        currentThick=(i-1+0.5)/fiberNum*totalThick;
        if currentThick<tempThick
            S=S+currentThick*totalThick/fiberNum*tempE;
            A=A+totalThick/fiberNum*tempE;
        else
            tempLayer=tempLayer+1;
            tempThick=tempThick+thickness(tempLayer);
            tempE=E(tempLayer);
            S=S+currentThick*totalThick/fiberNum*tempE;
            A=A+totalThick/fiberNum*tempE;
        end
    end
    
    % Solve for the centroid
    C=S/A;
    
    % Solve for the second moment of inertia
    EIcomp=0;
    tempThick=thickness(1);    
    tempE=E(1);
    tempLayer=1;
    
    for i=1:fiberNum
        currentThick=(i-1+0.5)/fiberNum*totalThick;
        if currentThick<tempThick
            EIcomp=EIcomp+(currentThick-C)^2*totalThick/fiberNum*tempE;
        else
            tempLayer=tempLayer+1;
            tempThick=tempThick+thickness(tempLayer);
            tempE=E(tempLayer);
            EIcomp=EIcomp+(currentThick-C)^2*totalThick/fiberNum*tempE;
        end
    end   
        

end


