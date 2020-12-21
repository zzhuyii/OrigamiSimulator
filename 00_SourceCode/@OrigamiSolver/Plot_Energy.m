%% The following code is used to plot the energy versus loading process
%
% Plot the strain energy history of the loading and self-assemble process.
% The function takes in input for both processes
%
% Input: 
%       UhisLoading: deformation history for loading;
%       strainEnergyLoading: strain energy history for loading;
%       UhisAssemble: deformaito history for self-assemble;
%       strainEnergyAssemble: strain energy history for assemble process;
%

function Plot_Energy(obj,UhisLoading,strainEnergyLoading)
    A=size(UhisLoading);
    IncreStep=A(1);
    UhisLoadingContract=zeros(IncreStep,1);
    Step=zeros(IncreStep,1);
    for i=1:IncreStep
        Step(i)=i;
        tempU=squeeze(UhisLoading(i,:,:));
        UhisLoadingContract(i)=norm(tempU) ;
        %UhisLoading(i)=norm(tempU(30,3));
    end

    % Plot Dots
    StrainEnergy=strainEnergyLoading;
    
    figure
    Strain1=StrainEnergy(:,1);
    Strain4=StrainEnergy(:,1)+StrainEnergy(:,4);
    Strain2=StrainEnergy(:,1)+StrainEnergy(:,2)+StrainEnergy(:,4);
    Strain3=StrainEnergy(:,1)+StrainEnergy(:,2)+StrainEnergy(:,3)+StrainEnergy(:,4);

    hold on
    
%     plot(UhisLoading,Strain1);
%     plot(UhisLoading,Strain4);
%     plot(UhisLoading,Strain2);
%     plot(UhisLoading,Strain3);

    plot(Step,Strain1);
    plot(Step,Strain4);
    plot(Step,Strain2);
    plot(Step,Strain3);

    legend({'Crease Bending','Crease Stretching','Panel Bending','Panel Stretching'},'Location','northwest');
    hold off
    
end