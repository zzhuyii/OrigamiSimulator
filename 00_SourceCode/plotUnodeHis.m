%% Plot the loading history

function plotUnodeHis(Uhis,LoadConstant)

IncreStep=LoadConstant(1);
h=LoadConstant(2);

A=size(Uhis);
IncreStep=A(1);
n1=A(2);
n2=A(3);
UhisLoading=zeros(IncreStep,1);
time=zeros(IncreStep,1);
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
    time(i)=i*h;
end

% Plot Dots
figure
plot(time,UhisLoading,'--o');xlabel('time');ylabel('Disp History of selected node');


