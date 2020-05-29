%% Generate geometry of Miura Tube
% This code is used to generate the geometry of miura tubes if needed

function [Node,Panel] = GenerateMiuraTube(a,b1,b2,gama1,gama2,N,sequence,Extention)

Node=zeros(4+8*N,3);

l11=sqrt((a-b1*cos(gama1))^2+(b1*sin(gama1))^2);
l12=sqrt((b1*cos(gama1)+a)^2+(b1*sin(gama1))^2);
l21=sqrt((a-b2*cos(gama2))^2+(b2*sin(gama2))^2);
l22=sqrt((b2*cos(gama2)+a)^2+(b2*sin(gama2))^2);

h1=b1*sin(gama1)*Extention;
theta1=asin(sin(gama1)*Extention);
temp=(l12^2-h1^2-a^2-(b1*cos(theta1))^2)/(2*a*b1*cos(theta1));
beta=acos(temp);
sqrt(h1^2+(a*sin(beta))^2+(a*cos(beta)+b1*cos(theta1))^2);

deltaZ=a*sin(beta);
theta2=acos((b2^2+(a*cos(beta))^2-(l21^2-deltaZ^2))/(2*b2*a*cos(beta)));

deltaX1=b1*sin(theta1);
deltaX2=b2*sin(theta2);

deltaY10=b1*cos(theta1);
deltaY11=a*cos(beta);
deltaY20=b2*cos(theta2);
deltaY21=a*cos(beta);

% Generate the Node Location
if sequence(1)==1 
    Node(1,:)=[0 0 0];
    Node(2,:)=[0 deltaY11 deltaZ];
    Node(3,:)=[0 2*deltaY11 0];
    Node(4,:)=[0 deltaY11 -deltaZ];
else
    Node(1,:)=[0 0 0];
    Node(2,:)=[0 deltaY21 deltaZ];
    Node(3,:)=[0 2*deltaY21 0];
    Node(4,:)=[0 deltaY21 -deltaZ];
end

X=0;

for i=1:N
    if sequence(i)==1
        X=X+deltaX1;        
        Node((8*i-4)+1,:)=[X deltaY10 0];
        Node((8*i-4)+2,:)=[X deltaY10+deltaY11 deltaZ];
        Node((8*i-4)+3,:)=[X deltaY10+2*deltaY11 0];
        Node((8*i-4)+4,:)=[X deltaY10+deltaY11 -deltaZ];
        X=X+deltaX1;
        Node((8*i-4)+5,:)=[X 0 0];
        Node((8*i-4)+6,:)=[X deltaY11 deltaZ];
        Node((8*i-4)+7,:)=[X 2*deltaY11 0];
        Node((8*i-4)+8,:)=[X deltaY11 -deltaZ];
    else
        X=X+deltaX2;        
        Node((8*i-4)+1,:)=[X deltaY20 0];
        Node((8*i-4)+2,:)=[X deltaY20+deltaY21 deltaZ];
        Node((8*i-4)+3,:)=[X deltaY20+2*deltaY21 0];
        Node((8*i-4)+4,:)=[X deltaY20+deltaY21 -deltaZ];
        X=X+deltaX2;
        Node((8*i-4)+5,:)=[X 0 0];
        Node((8*i-4)+6,:)=[X deltaY21 deltaZ];
        Node((8*i-4)+7,:)=[X 2*deltaY21 0];
        Node((8*i-4)+8,:)=[X deltaY21 -deltaZ];
    end
end

% Generate Panels
for i=1:N
   base=8*(i-1)
   Panel{base+1}=[base+1 base+5 base+6 base+2];
   Panel{base+2}=[base+2 base+6 base+7 base+3];
   Panel{base+3}=[base+3 base+7 base+8 base+4];
   Panel{base+4}=[base+4 base+8 base+5 base+1];
   Panel{base+5}=[base+5 base+9 base+10 base+6];
   Panel{base+6}=[base+6 base+10 base+11 base+7];
   Panel{base+7}=[base+7 base+11 base+12 base+8];
   Panel{base+8}=[base+8 base+12 base+9 base+5];
end




