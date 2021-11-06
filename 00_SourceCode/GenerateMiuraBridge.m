%% Generate geometry of Miura sheet
% This code is used to generate the geometry of miura sheet if needed

function[Node,Panel] =GenerateMiuraBridge(a,b,gama,m,n,Ext)
h=b*sin(gama)*Ext;
theta=asin(Ext*sin(gama));
l=sqrt((a+b*cos(gama))^2+(b*sin(gama))^2);
temp=(l^2-h^2-a^2-(b*cos(theta))^2)/(2*a*b*cos(theta));
beta=norm(acos(temp));
Beta=beta*180/3.14

deltaZ=a*sin(beta);
deltaX=b*sin(theta);
deltaY0=b*cos(theta);
deltaY1=a*cos(beta);

panelNum=m*n;
nodeNum=(m+1)*(n+1);
Panel{1}=zeros(1);

for i=1:m+1
    for j=1:n+1
        num=(i-1)*(n+1)+j;
        initial=deltaY0;
        if mod(i,2)
           initial=0; 
        end
        if mod(j,2)
            Node(num,:)=[(i-1)*deltaX (j-1)*deltaY1+initial 0];
        else
            Node(num,:)=[(i-1)*deltaX (j-1)*deltaY1+initial deltaZ]; 
        end       
    end    
end

for i=1:m
    for j=1:n
       Panel{(i-1)*n+j}=[(i-1)*(n+1)+j i*(n+1)+j i*(n+1)+j+1 (i-1)*(n+1)+j+1];
    end
end
    

