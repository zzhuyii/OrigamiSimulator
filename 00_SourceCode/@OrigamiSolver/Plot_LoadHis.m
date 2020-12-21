%% Plot the loading history
%
% The function plot the force-displacement curve. Various displacement
% quantities can be used for the plotting. Please double check the setup.
%
% Input:
%       loadHis: loading history;
%       UhisLoading: deformation history of the loading process;
%

function Plot_LoadHis(obj,loadHis,UhisLoading)

A=size(UhisLoading);
IncreStep=A(1);
UhisLoadingContract=zeros(IncreStep,1);
for i=1:IncreStep
    tempU=squeeze(UhisLoading(i,:,:));
    %tempU1=zeros(1,3);
    %tempU1(1,:)=Uhis(1,46,:);
    UhisLoadingContract(i)=norm(tempU);    
    %UhisLoading(i)=tempU(6,3);
end
% Plot Dots
figure
plot(UhisLoadingContract,loadHis,'--o');xlabel('Displacement History');ylabel('Load History');


