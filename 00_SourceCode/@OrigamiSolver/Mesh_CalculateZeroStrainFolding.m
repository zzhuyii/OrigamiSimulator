%% Calculate the zero strain folding angle of rotational springs

function sprZeroStrain=Mesh_CalculateZeroStrainFolding(obj,rotZeroStrain)
  
    A=size(obj.barLength);
    barNum=A(1);
    sprZeroStrain=pi*ones(barNum,1);    
    
    % when we need to consider the compliant creases
    if obj.compliantCreaseOpen==1
        % calculate the properties for bar/spr in the crease region
        for i=1:obj.oldCreaseNum
           if obj.oldCreaseType(i)>1
               theta0=rotZeroStrain(i); 
               
               % asign targe zero strain of bars
               sprZeroStrain(obj.creaseRef(i,1))=...
                   (theta0-pi)*(1-obj.zeroStrianDistributionFactor)/2+pi;
               sprZeroStrain(obj.creaseRef(i,2))=...
                   (theta0-pi)*(1-obj.zeroStrianDistributionFactor)/2+pi;
               sprZeroStrain(obj.creaseRef(i,5))=...
                   (theta0-pi)*(obj.zeroStrianDistributionFactor)+pi;
               sprZeroStrain(obj.creaseRef(i,6))=...
                   (theta0-pi)*(obj.zeroStrianDistributionFactor)+pi;               
            
           end
        end
        
    % when we donot need to consider the compliant creases
    else
        % Assemble the stiffness properties for creases
        for i=1:obj.oldCreaseNum
            if obj.oldCreaseType(i)>1
                theta0=rotZeroStrain(i);     
                sprZeroStrain(i)=theta0;
            end
        end
    end
end
