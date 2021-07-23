%% Generate geometry of Yoshiura sheet
% This code is used to generate the geometry of miura sheet if needed
% both n and m should be even number

function[Node,Panel] =GenerateYoshimuraSheet(a,b,m,n,Ext)

delta=b*sqrt((1-Ext^2));
sai=atan(a/2/delta);
theta=pi-2*sai;
L=a/2/sin(theta);

nodeNum=m*(2*n+1)/2+n;
Node=zeros(nodeNum,3);
Panel{1}=zeros(1);

for i=1:(m/2)    
    for j=1:n
        Node((i-1)*(2*n+1)+j,:)=...
            [b*Ext*(2*(i-1)), ...
             L*sin((2*j-n-1)*theta),L*cos((2*j-n-1)*theta)];
    end        
    for j=1:n+1
        Node((i-1)*(2*n+1)+j+n,:)=...
            [b*Ext*(2*(i-1)+1), ...
             L*sin((2*j-n-2)*theta),L*cos((2*j-n-2)*theta)];
    end    
end

for j=1:n
    Node(m/2*(2*n+1)+j,:)=...
        [b*Ext*(m), ...
         L*sin((2*j-n-1)*theta),L*cos((2*j-n-1)*theta)];
end    

Node(:,3)=Node(:,3)-max(Node(:,3));

for i=1:(m/2)
    for j=1:n-1
       Panel{(i-1)*(4*n-2)+j}=[(i-1)*(2*n+1)+j ...
           (i-1)*(2*n+1)+n+1+j  (i-1)*(2*n+1)+j+1];
    end
    for j=1:n
       Panel{(i-1)*(4*n-2)+n-1+j}=[(i-1)*(2*n+1)+n+j ...
           (i-1)*(2*n+1)+j  (i-1)*(2*n+1)+n+j+1];
    end
    for j=1:n
       Panel{(i-1)*(4*n-2)+2*n-1+j}=[(i-1)*(2*n+1)+n+j ...
           (i)*(2*n+1)+j  (i-1)*(2*n+1)+n+j+1];
    end
    for j=1:n-1
       Panel{(i-1)*(4*n-2)+3*n-1+j}=[(i)*(2*n+1)+j ...
           (i-1)*(2*n+1)+n+1+j  (i)*(2*n+1)+j+1];
    end
end
    



