%% This function add a rotational spring for building thick panels
% The funciton need to specify 8 nodes to work. These 8 nodes are from the
% adjacent panels where the hinge will be added. 

function AddHingeForThickPanel(obj,n1,n2,n3,n4,n5,n6,n7,n8,BarArea,BarAreaLong,L,Gap,t,hingeStiff)


    A=size(obj.newNode);
    nodeNum=A(1);
    A=size(obj.barType);
    barNum=A(1);    
    
    obj.newNode(nodeNum+1,:)=(obj.newNode(n1,:)+obj.newNode(n2,:))/2;
    obj.newNode(nodeNum+2,:)=(obj.newNode(n3,:)+obj.newNode(n4,:))/2;
    
    obj.barType(barNum+1)=1;
    obj.barConnect(barNum+1,:)=[n1,nodeNum+1];
    obj.barArea(barNum+1)=BarArea;
    obj.barLength(barNum+1)=Gap/2;        
    
    obj.barType(barNum+2)=1;
    obj.barConnect(barNum+2,:)=[n3,nodeNum+2];
    obj.barArea(barNum+2)=BarArea;
    obj.barLength(barNum+2)=Gap/2;        
    
    obj.barType(barNum+3)=1;
    obj.barConnect(barNum+3,:)=[n2,nodeNum+1];
    obj.barArea(barNum+3)=BarArea;
    obj.barLength(barNum+3)=Gap/2;        
    
    obj.barType(barNum+4)=1;
    obj.barConnect(barNum+4,:)=[n4,nodeNum+2];
    obj.barArea(barNum+4)=BarArea;
    obj.barLength(barNum+4)=Gap/2;      
    
    obj.barType(barNum+5)=1;
    obj.barConnect(barNum+5,:)=[n5,nodeNum+1];
    obj.barArea(barNum+5)=BarArea;
    obj.barLength(barNum+5)=sqrt((Gap/2)^2+t^2);       
    
    obj.barType(barNum+6)=1;
    obj.barConnect(barNum+6,:)=[n7,nodeNum+2];
    obj.barArea(barNum+6)=BarArea;
    obj.barLength(barNum+6)=sqrt((Gap/2)^2+t^2);     
    
    obj.barType(barNum+7)=1;
    obj.barConnect(barNum+7,:)=[n6,nodeNum+1];
    obj.barArea(barNum+7)=BarArea;
    obj.barLength(barNum+7)=sqrt((Gap/2)^2+t^2);       
    
    obj.barType(barNum+8)=1;
    obj.barConnect(barNum+8,:)=[n8,nodeNum+2];
    obj.barArea(barNum+8)=BarArea;
    obj.barLength(barNum+8)=sqrt((Gap/2)^2+t^2);     
    
    
    %% Cross bars
    obj.barType(barNum+9)=1;
    obj.barConnect(barNum+9,:)=[n5,nodeNum+2];
    obj.barArea(barNum+9)=BarAreaLong;
    obj.barLength(barNum+9)=sqrt(L^2+t*t+Gap*Gap/4);       
    
    obj.barType(barNum+10)=1;
    obj.barConnect(barNum+10,:)=[n7,nodeNum+1];
    obj.barArea(barNum+10)=BarAreaLong;
    obj.barLength(barNum+10)=sqrt(L^2+t*t+Gap*Gap/4);    
    
    obj.barType(barNum+11)=1;
    obj.barConnect(barNum+11,:)=[n6,nodeNum+2];
    obj.barArea(barNum+11)=BarAreaLong;
    obj.barLength(barNum+11)=sqrt(L^2+t*t+Gap*Gap/4);          
    
    obj.barType(barNum+12)=1;
    obj.barConnect(barNum+12,:)=[n8,nodeNum+1];
    obj.barArea(barNum+12)=BarAreaLong;
    obj.barLength(barNum+12)=sqrt(L^2+t*t+Gap*Gap/4);     
    
    obj.barType(barNum+13)=1;
    obj.barConnect(barNum+13,:)=[n1,nodeNum+2];
    obj.barArea(barNum+13)=BarAreaLong;
    obj.barLength(barNum+13)=sqrt(L^2+Gap*Gap/4);         
    
    obj.barType(barNum+14)=1;
    obj.barConnect(barNum+14,:)=[n3,nodeNum+1];
    obj.barArea(barNum+14)=BarAreaLong;
    obj.barLength(barNum+14)=sqrt(L^2+Gap*Gap/4);     
    
    obj.barType(barNum+15)=1;
    obj.barConnect(barNum+15,:)=[n2,nodeNum+2];
    obj.barArea(barNum+15)=BarAreaLong;
    obj.barLength(barNum+15)=sqrt(L^2+Gap*Gap/4);         
    
    obj.barType(barNum+16)=1;
    obj.barConnect(barNum+16,:)=[n4,nodeNum+1];
    obj.barArea(barNum+16)=BarAreaLong;
    obj.barLength(barNum+16)=sqrt(L^2+Gap*Gap/4);    
    
    % Final rotational hinge
    obj.barType(barNum+17)=1;
    obj.barConnect(barNum+17,:)=[nodeNum+1,nodeNum+2];
    obj.barArea(barNum+17)=BarAreaLong;
    obj.barLength(barNum+17)=L; 
    obj.sprK(barNum+17)=hingeStiff;
    obj.sprIJKL(barNum+17,:)=[n1,nodeNum+1,nodeNum+2,n4]; 

end