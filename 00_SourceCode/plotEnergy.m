%% The following code is used to plot the energy versus loading process

function plotEnergy(Uhis,StrainEnergyLoading,assembleHis,StrainEnergyAssemble)
    A=size(Uhis);
    IncreStep=A(1);
    n1=A(2);
    n2=A(3);
    UhisLoading=zeros(IncreStep,1);
    for i=1:IncreStep
        tempU=zeros(n1,n2);
        for j=1:n1
            for k=1:n2
               tempU(j,k)=Uhis(i,j,k);
            end
        end
        UhisLoading(i)=norm(tempU) ;
        %UhisLoading(i)=norm(tempU(30,3));
    end

    A=size(assembleHis);
    IncreStep=A(1);
    n1=A(2);
    n2=A(3);
    UhisAssemble=zeros(IncreStep,1);
    for i=1:IncreStep
        tempU=zeros(n1,n2);
        for j=1:n1
            for k=1:n2
               tempU(j,k)=assembleHis(i,j,k);
            end
        end
        UhisAssemble(i)=norm(tempU) ;
        %UhisAssemble(i)=tempU(30,3);
    end

    % Plot Dots
    Uhis=[UhisAssemble;UhisLoading;];
    N=size(Uhis,1);
    Step=zeros(N,1);
    for i=1:N
        Step(i)=i;
    end

    StrainEnergy=[StrainEnergyAssemble;StrainEnergyLoading;];

    figure
    Strain1=StrainEnergy(:,1);
    Strain4=StrainEnergy(:,1)+StrainEnergy(:,4);
    Strain2=StrainEnergy(:,1)+StrainEnergy(:,2)+StrainEnergy(:,4);
    Strain3=StrainEnergy(:,1)+StrainEnergy(:,2)+StrainEnergy(:,3)+StrainEnergy(:,4);

    hold on
    plot(Uhis,Strain1);
    plot(Uhis,Strain4);
    plot(Uhis,Strain2);
    plot(Uhis,Strain3);

%     plot(Step,Strain1);
%     plot(Step,Strain4);
%     plot(Step,Strain2);
%     plot(Step,Strain3);
    legend({'Crease Bending','Crease Stretching','Panel Bending','Panel Stretching'},'Location','northwest');
    hold off
    
end