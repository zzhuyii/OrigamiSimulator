%% This code generates the thick Yoshimura Pattern

function ori=GenerateThickYoshimuraUnitCell(L,t,gap)

ori=OrigamiSolver;

ori.panelInnerBarStart=1;
panelBarArea=100;
connectorBarArea=100;
hingeStiff=1000;

%% Add nodes bottom layer
% Add node for the first row of panels
ori.newNode(1,:)=[(sqrt(2)/2)*gap*1+L*1,L,0];
ori.newNode(2,:)=[(sqrt(2)/2)*gap*1+L*1,0,0];
ori.newNode(3,:)=[(sqrt(2)/2)*gap*1+L*0,0,0];

ori.newNode(4,:)=[0,(sqrt(2)/2)*gap,0];
ori.newNode(5,:)=[0,(sqrt(2)/2)*gap+gap+2*L,0];
ori.newNode(6,:)=[0,(sqrt(2)/2)*gap+gap/2+L,0];
ori.newNode(7,:)=[L,L+(sqrt(2)/2)*gap,0];
ori.newNode(8,:)=[L,L+(sqrt(2)/2)*gap+gap,0];

ori.newNode(9,:)=[(sqrt(2)/2)*gap,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(10,:)=[(sqrt(2)/2)*gap+L,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(11,:)=[(sqrt(2)/2)*gap+L,(sqrt(2)/2)*gap*2+gap+L*1,0]; 


% Add node for the second row of panels
ori.newNode(12,:)=[(sqrt(2)/2)*gap+L+gap,L,0];
ori.newNode(13,:)=[(sqrt(2)/2)*gap+L+gap,0,0];
ori.newNode(14,:)=[(sqrt(2)/2)*gap+gap+L*2,0,0];

ori.newNode(15,:)=[2*L+sqrt(2)*gap+gap,(sqrt(2)/2)*gap,0];
ori.newNode(16,:)=[2*L+sqrt(2)*gap+gap,(sqrt(2)/2)*gap+gap+2*L,0];
ori.newNode(17,:)=[L*2+gap+sqrt(2)*gap,(sqrt(2)/2)*gap+gap/2+L,0];
ori.newNode(18,:)=[L+sqrt(2)*gap+gap,L+(sqrt(2)/2)*gap,0];
ori.newNode(19,:)=[L+sqrt(2)*gap+gap,L+(sqrt(2)/2)*gap+gap,0];

ori.newNode(20,:)=[(sqrt(2)/2)*gap+gap+2*L,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(21,:)=[(sqrt(2)/2)*gap+L+gap,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(22,:)=[(sqrt(2)/2)*gap+L+gap,(sqrt(2)/2)*gap*2+gap+L*1,0]; 

%% Add nodes top layer
% Add node for the first row of panels
ori.newNode(22+1,:)=[(sqrt(2)/2)*gap*1+L*1,L,0];
ori.newNode(22+2,:)=[(sqrt(2)/2)*gap*1+L*1,0,0];
ori.newNode(22+3,:)=[(sqrt(2)/2)*gap*1+L*0,0,0];

ori.newNode(22+4,:)=[0,(sqrt(2)/2)*gap,0];
ori.newNode(22+5,:)=[0,(sqrt(2)/2)*gap+gap+2*L,0];
ori.newNode(22+6,:)=[0,(sqrt(2)/2)*gap+gap/2+L,0];
ori.newNode(22+7,:)=[L,L+(sqrt(2)/2)*gap,0];
ori.newNode(22+8,:)=[L,L+(sqrt(2)/2)*gap+gap,0];

ori.newNode(22+9,:)=[(sqrt(2)/2)*gap,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(22+10,:)=[(sqrt(2)/2)*gap+L,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(22+11,:)=[(sqrt(2)/2)*gap+L,(sqrt(2)/2)*gap*2+gap+L*1,0]; 


% Add node for the second row of panels
ori.newNode(22+12,:)=[(sqrt(2)/2)*gap+L+gap,L,0];
ori.newNode(22+13,:)=[(sqrt(2)/2)*gap+L+gap,0,0];
ori.newNode(22+14,:)=[(sqrt(2)/2)*gap+gap+L*2,0,0];

ori.newNode(22+15,:)=[2*L+sqrt(2)*gap+gap,(sqrt(2)/2)*gap,0];
ori.newNode(22+16,:)=[2*L+sqrt(2)*gap+gap,(sqrt(2)/2)*gap+gap+2*L,0];
ori.newNode(22+17,:)=[L*2+gap+sqrt(2)*gap,(sqrt(2)/2)*gap+gap/2+L,0];
ori.newNode(22+18,:)=[L+sqrt(2)*gap+gap,L+(sqrt(2)/2)*gap,0];
ori.newNode(22+19,:)=[L+sqrt(2)*gap+gap,L+(sqrt(2)/2)*gap+gap,0];

ori.newNode(22+20,:)=[(sqrt(2)/2)*gap+gap+2*L,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(22+21,:)=[(sqrt(2)/2)*gap+L+gap,(sqrt(2)/2)*gap*2+gap+L*2,0]; 
ori.newNode(22+22,:)=[(sqrt(2)/2)*gap+L+gap,(sqrt(2)/2)*gap*2+gap+L*1,0]; 

ori.newNode(23:44,3)=ori.newNode(23:44,3)+t;


%% Add Panels for Plotting
% Panel of the faces
ori.newPanel{1}=[1,2,3];
ori.newPanel{2}=[4,5,8,7];
ori.newPanel{3}=[9,10,11];

ori.newPanel{4}=[12,13,14];
ori.newPanel{5}=[15,16,19,18];
ori.newPanel{6}=[20,21,22];

ori.newPanel{7}=[23,24,25];
ori.newPanel{8}=[26,27,30,29];
ori.newPanel{9}=[31,32,33];

ori.newPanel{10}=[34,35,36];
ori.newPanel{11}=[37,38,41,40];
ori.newPanel{12}=[42,43,44];

% panel for the side faces
ori.newPanel{13}=[1,2,2+22,1+22];
ori.newPanel{14}=[2,3,3+22,2+22];
ori.newPanel{15}=[1,3,3+22,1+22];

ori.newPanel{16}=[4,5,5+22,4+22];
ori.newPanel{17}=[5,8,8+22,5+22];
ori.newPanel{18}=[8,7,7+22,8+22];
ori.newPanel{19}=[4,7,7+22,4+22];

ori.newPanel{20}=[9,10,10+22,9+22];
ori.newPanel{21}=[10,11,11+22,10+22];
ori.newPanel{22}=[9,11,11+22,9+22];

% Second row
ori.newPanel{23}=[1+11,2+11,2+22+11,1+22+11];
ori.newPanel{24}=[2+11,3+11,3+22+11,2+22+11];
ori.newPanel{25}=[1+11,3+11,3+22+11,1+22+11];

ori.newPanel{26}=[4+11,5+11,5+22+11,4+22+11];
ori.newPanel{27}=[5+11,8+11,8+22+11,5+22+11];
ori.newPanel{28}=[8+11,7+11,7+22+11,8+22+11];
ori.newPanel{29}=[4+11,7+11,7+22+11,4+22+11];

ori.newPanel{30}=[9+11,10+11,10+22+11,9+22+11];
ori.newPanel{31}=[10+11,11+11,11+22+11,10+22+11];
ori.newPanel{32}=[9+11,11+11,11+22+11,9+22+11];

%% Add bars for bottom layer
% Add bars for the first row of panels
ori.AddBar(1,2,panelBarArea,L);
ori.AddBar(2,3,panelBarArea,L);
ori.AddBar(1,3,panelBarArea,L*sqrt(2));

ori.AddBar(4,6,panelBarArea,L+gap/2);
ori.AddBar(5,6,panelBarArea,L+gap/2);
ori.AddBar(4,7,panelBarArea,L*sqrt(2));
ori.AddBar(5,8,panelBarArea,L*sqrt(2));

ori.AddBar(7,8,panelBarArea,gap);
ori.AddBar(6,7,panelBarArea,sqrt(L^2+(gap/2)^2));
ori.AddBar(6,8,panelBarArea,sqrt(L^2+(gap/2)^2));

ori.AddBar(9,10,panelBarArea,L);
ori.AddBar(10,11,panelBarArea,L);
ori.AddBar(9,11,panelBarArea,L*sqrt(2));

% Add bars for the second row of panels
ori.AddBar(12,13,panelBarArea,L);
ori.AddBar(13,14,panelBarArea,L);
ori.AddBar(12,14,panelBarArea,L*sqrt(2));

ori.AddBar(15,17,panelBarArea,L+gap/2);
ori.AddBar(16,17,panelBarArea,L+gap/2);
ori.AddBar(15,18,panelBarArea,L*sqrt(2));
ori.AddBar(16,19,panelBarArea,L*sqrt(2));

ori.AddBar(18,19,panelBarArea,gap);
ori.AddBar(17,18,panelBarArea,sqrt(L^2+(gap/2)^2));
ori.AddBar(17,19,panelBarArea,sqrt(L^2+(gap/2)^2));

ori.AddBar(22,20,panelBarArea,L);
ori.AddBar(20,21,panelBarArea,L);
ori.AddBar(22,21,panelBarArea,L*sqrt(2));


%% Add bars for top layer
% Add bars for the first row of panels
ori.AddBar(22+1,22+2,panelBarArea,L);
ori.AddBar(22+2,22+3,panelBarArea,L);
ori.AddBar(22+1,22+3,panelBarArea,L*sqrt(2));

ori.AddBar(22+4,22+6,panelBarArea,L+gap/2);
ori.AddBar(22+5,22+6,panelBarArea,L+gap/2);
ori.AddBar(22+4,22+7,panelBarArea,L*sqrt(2));
ori.AddBar(22+5,22+8,panelBarArea,L*sqrt(2));

ori.AddBar(22+7,22+8,panelBarArea,gap);
ori.AddBar(22+6,22+7,panelBarArea,sqrt(L^2+(gap/2)^2));
ori.AddBar(22+6,22+8,panelBarArea,sqrt(L^2+(gap/2)^2));

ori.AddBar(22+9,22+10,panelBarArea,L);
ori.AddBar(22+10,22+11,panelBarArea,L);
ori.AddBar(22+9,22+11,panelBarArea,L*sqrt(2));

% Add bars for the second row of panels
ori.AddBar(22+12,22+13,panelBarArea,L);
ori.AddBar(22+13,22+14,panelBarArea,L);
ori.AddBar(22+12,22+14,panelBarArea,L*sqrt(2));

ori.AddBar(22+15,22+17,panelBarArea,L+gap/2);
ori.AddBar(22+16,22+17,panelBarArea,L+gap/2);
ori.AddBar(22+15,22+18,panelBarArea,L*sqrt(2));
ori.AddBar(22+16,22+19,panelBarArea,L*sqrt(2));

ori.AddBar(22+18,22+19,panelBarArea,gap);
ori.AddBar(22+17,22+18,panelBarArea,sqrt(L^2+(gap/2)^2));
ori.AddBar(22+17,22+19,panelBarArea,sqrt(L^2+(gap/2)^2));

ori.AddBar(22+22,22+20,panelBarArea,L);
ori.AddBar(22+20,22+21,panelBarArea,L);
ori.AddBar(22+22,22+21,panelBarArea,L*sqrt(2));

%% Add vertical connection bars
for i=1:22
    ori.AddBar(i,22+i,panelBarArea,t);
end


%% Add diagonal bars in vertical direction
% Add bars for the first row of panels
ori.AddBar(1,22+2,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(2,22+3,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(1,22+3,panelBarArea,sqrt(4*L^2+t^2));

ori.AddBar(4,22+6,panelBarArea,sqrt((L+gap/2)^2+t^2));
ori.AddBar(5,22+6,panelBarArea,sqrt((L+gap/2)^2+t^2));
ori.AddBar(4,22+7,panelBarArea,sqrt(4*L^2+t^2));
ori.AddBar(5,22+8,panelBarArea,sqrt(4*L^2+t^2));

ori.AddBar(7,22+8,panelBarArea,sqrt(gap^2+t^2));
ori.AddBar(6,22+7,panelBarArea,sqrt(L^2+(gap/2)^2+t^2));
ori.AddBar(6,22+8,panelBarArea,sqrt(L^2+(gap/2)^2+t^2));

ori.AddBar(9,22+10,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(10,22+11,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(9,22+11,panelBarArea,sqrt(4*L^2+t^2));

% Add bars for the second row of panels
ori.AddBar(12,22+13,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(13,22+14,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(12,22+14,panelBarArea,sqrt(4*L^2+t^2));

ori.AddBar(15,22+17,panelBarArea,sqrt((L+gap/2)^2+t^2));
ori.AddBar(16,22+17,panelBarArea,sqrt((L+gap/2)^2+t^2));
ori.AddBar(15,22+18,panelBarArea,sqrt(4*L^2+t^2));
ori.AddBar(16,22+19,panelBarArea,sqrt(4*L^2+t^2));

ori.AddBar(18,22+19,panelBarArea,sqrt(gap^2+t^2));
ori.AddBar(17,22+18,panelBarArea,sqrt(L^2+(gap/2)^2+t^2));
ori.AddBar(17,22+19,panelBarArea,sqrt(L^2+(gap/2)^2+t^2));

ori.AddBar(22,22+20,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(20,22+21,panelBarArea,sqrt(L^2+t^2));
ori.AddBar(22,22+21,panelBarArea,sqrt(4*L^2+t^2));


%% Add hinge for connecting thick panel    
creaseNumBeforeHinge=length(ori.sprK);

ori.AddHingeForThickPanel(23,29,25,26,1,7,3,4,connectorBarArea,connectorBarArea/10,sqrt(2)*L,gap,t,hingeStiff);
ori.AddHingeForThickPanel(38,42,41,44,16,20,19,22,connectorBarArea,connectorBarArea/10,sqrt(2)*L,gap,t,hingeStiff);

ori.AddHingeForThickPanel(30,33,27,31,8,11,5,9,connectorBarArea,connectorBarArea/10,sqrt(2)*L,gap,t,hingeStiff);
ori.AddHingeForThickPanel(36,37,34,40,14,15,12,18,connectorBarArea,connectorBarArea/10,sqrt(2)*L,gap,t,hingeStiff);

ori.AddHingeForThickPanel(2,13,1,12,24,35,23,34,connectorBarArea,connectorBarArea/10,L,gap,t,hingeStiff);
ori.AddHingeForThickPanel(11,22,10,21,33,44,32,43,connectorBarArea,connectorBarArea/10,L,gap,t,hingeStiff);

creaseNumAfterHinge=length(ori.sprK); 

%% Set up terms for initialization
newNodeNum = size(ori.newNode);
newNodeNum = newNodeNum(1);
ori.currentAppliedForce = zeros(newNodeNum,3);   
ori.currentU = zeros(newNodeNum,3);
ori.currentSprZeroStrain = pi*logical(ori.sprK);
ori.currentRotZeroStrain = pi*logical(ori.sprK);

ori.oldCreaseNum=length(squeeze(ori.sprK));
ori.oldCreaseType=zeros(creaseNumAfterHinge,1);
ori.oldCreaseType(creaseNumBeforeHinge:creaseNumAfterHinge)=2;


end

