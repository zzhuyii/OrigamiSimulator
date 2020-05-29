%% Plot the loading history

function plotLoadHis(loadHis,Uhis)

A=size(Uhis);
IncreStep=A(1);
n1=A(2);
n2=A(3);
UhisLoading=zeros(IncreStep,1);
for i=1:IncreStep
    tempU=zeros(n1,n2);
    %tempU1=zeros(1,3);
    %tempU1(1,:)=Uhis(1,46,:);
    for j=1:n1
        for k=1:n2
           tempU(j,k)=Uhis(i,j,k);
        end
    end
    UhisLoading(i)=norm(tempU);    
    UhisLoading(i)=tempU(6,3);
end
% Plot Dots
figure
plot(UhisLoading,loadHis,'--o');xlabel('Displacement History');ylabel('Load History');


