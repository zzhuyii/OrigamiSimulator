%% This code generates the thick Yoshimura Pattern

function ori=GenerateModularOrigami(L,N,t,f,gap)

ori=OrigamiSolver;

ori.panelInnerBarStart=1;
panelBarArea=0.000258/3;
connectorBarArea=0.000258;

hingeStiff=0.000001;
hingeStiffYoshimura=f*hingeStiff;



% Initialize the size
ori.barArea=zeros(2,1);
ori.barLength=zeros(2,1);
ori.barType=zeros(2,1);
ori.sprK=zeros(2,1);

for i=1:N    

    for j =1:4        
        
        
        %% First Plate        
        if (mod(i,2)==1 && mod(j,2)==1) || (mod(i,2)==0 && mod(j,2)==0)
            ori.newNode(48*(i-1)+12*(j-1)+1,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1)+0,gap*i*(1+sqrt(2)/2)+L*(i-1)+0,0];
            ori.newNode(48*(i-1)+12*(j-1)+2,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1),gap*i*(1+sqrt(2)/2)+L*(i-1)+L,0];
            ori.newNode(48*(i-1)+12*(j-1)+3,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1)+L,gap*i*(1+sqrt(2)/2)+L*(i-1)+L,0];
            ori.newNode(48*(i-1)+12*(j-1)+4,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1),gap*i*(1+sqrt(2)/2)+L*(i-1)+0,t];
            ori.newNode(48*(i-1)+12*(j-1)+5,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1),gap*i*(1+sqrt(2)/2)+L*(i-1)+L,t];
            ori.newNode(48*(i-1)+12*(j-1)+6,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1)+L,gap*i*(1+sqrt(2)/2)+L*(i-1)+L,t];
        else
            ori.newNode(48*(i-1)+12*(j-1)+1,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1)+0,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+0,0];
            ori.newNode(48*(i-1)+12*(j-1)+2,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1),gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+L,0];
            ori.newNode(48*(i-1)+12*(j-1)+3,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1)+L,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1),0];
            ori.newNode(48*(i-1)+12*(j-1)+4,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1),gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+0,t];
            ori.newNode(48*(i-1)+12*(j-1)+5,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1),gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+L,t];
            ori.newNode(48*(i-1)+12*(j-1)+6,:)=[(1+sqrt(2)/2)*gap*(j-1)+L*(j-1)+L,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1),t];
        end 

        %% Second Plate    
        if (mod(i,2)==1 && mod(j,2)==1) || (mod(i,2)==0 && mod(j,2)==0)
            ori.newNode(48*(i-1)+12*(j-1)+7,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+0,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+0,0]; 
            ori.newNode(48*(i-1)+12*(j-1)+8,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1),0];
            ori.newNode(48*(i-1)+12*(j-1)+9,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+L,0];
            ori.newNode(48*(i-1)+12*(j-1)+10,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1),gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+0,t]; 
            ori.newNode(48*(i-1)+12*(j-1)+11,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1),t];
            ori.newNode(48*(i-1)+12*(j-1)+12,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)-gap*sqrt(2)/2+L*(i-1)+L,t];
        else
            ori.newNode(48*(i-1)+12*(j-1)+7,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+0,gap*i*(1+sqrt(2)/2)+L*(i-1)+L,0]; 
            ori.newNode(48*(i-1)+12*(j-1)+8,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)+L*(i-1),0];
            ori.newNode(48*(i-1)+12*(j-1)+9,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)+L*(i-1)+L,0];
            ori.newNode(48*(i-1)+12*(j-1)+10,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1),gap*i*(1+sqrt(2)/2)+L*(i-1)+L,t]; 
            ori.newNode(48*(i-1)+12*(j-1)+11,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)+L*(i-1),t];
            ori.newNode(48*(i-1)+12*(j-1)+12,:)=[(1+sqrt(2)/2)*gap*(j-1)+gap*sqrt(2)/2+L*(j-1)+L,gap*i*(1+sqrt(2)/2)+L*(i-1)+L,t];
        end
        
        
        %% Assign Bar Info

        % for the bridge design, slab panels should have more area because
        % of the distributed plate. We assume that the addtional area is
        % three times the area of the side trusses
        if j==2
            panelBarAreaTemp=3*panelBarArea;
        elseif j==3
            panelBarAreaTemp=3*panelBarArea;
        else
            panelBarAreaTemp=panelBarArea;
        end


        for q=1:2
            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+1)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+1,:)=[48*(i-1)+12*(j-1)+6*(q-1)+1,48*(i-1)+12*(j-1)+6*(q-1)+2];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+1)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+1)=L;

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+2)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+2,:)=[48*(i-1)+12*(j-1)+6*(q-1)+2,48*(i-1)+12*(j-1)+6*(q-1)+3];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+2)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+2)=L;

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+3)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+3,:)=[48*(i-1)+12*(j-1)+6*(q-1)+1,48*(i-1)+12*(j-1)+6*(q-1)+3];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+3)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+3)=L*sqrt(2);

            % top triangle
            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+4)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+4,:)=[48*(i-1)+12*(j-1)+6*(q-1)+4,48*(i-1)+12*(j-1)+6*(q-1)+5];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+4)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+4)=L;

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+5)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+5,:)=[48*(i-1)+12*(j-1)+6*(q-1)+5,48*(i-1)+12*(j-1)+6*(q-1)+6];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+5)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+5)=L;

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+6)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+6,:)=[48*(i-1)+12*(j-1)+6*(q-1)+4,48*(i-1)+12*(j-1)+6*(q-1)+6];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+6)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+6)=L*sqrt(2);    

            % Vertical
            % Vertical bar represents the stiffness in the thickness
            % direction, because the area is much larger, we assume that
            % the vertical bar is 10 times the area
            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+7)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+7,:)=[48*(i-1)+12*(j-1)+6*(q-1)+1,48*(i-1)+12*(j-1)+6*(q-1)+4];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+7)=panelBarAreaTemp*10;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+7)=t;

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+8)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+8,:)=[48*(i-1)+12*(j-1)+6*(q-1)+2,48*(i-1)+12*(j-1)+6*(q-1)+5];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+8)=panelBarAreaTemp*10;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+8)=t;

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+9)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+9,:)=[48*(i-1)+12*(j-1)+6*(q-1)+3,48*(i-1)+12*(j-1)+6*(q-1)+6];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+9)=panelBarAreaTemp*10;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+9)=t;


            % Diagonal
            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+10)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+10,:)=[48*(i-1)+12*(j-1)+6*(q-1)+3,48*(i-1)+12*(j-1)+6*(q-1)+5];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+10)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+10)=sqrt(t*t+L*L);

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+11)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+11,:)=[48*(i-1)+12*(j-1)+6*(q-1)+2,48*(i-1)+12*(j-1)+6*(q-1)+4];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+11)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+11)=sqrt(t*t+L*L);

            ori.barType(96*(i-1)+24*(j-1)+12*(q-1)+12)=1;
            ori.barConnect(96*(i-1)+24*(j-1)+12*(q-1)+12,:)=[48*(i-1)+12*(j-1)+6*(q-1)+1,48*(i-1)+12*(j-1)+6*(q-1)+6];
            ori.barArea(96*(i-1)+24*(j-1)+12*(q-1)+12)=panelBarAreaTemp;
            ori.barLength(96*(i-1)+24*(j-1)+12*(q-1)+12)=sqrt(t*t+2*L*L);
        end

        %% Assign Triangle Info
        for q=1:2
            ori.newPanel{40*(i-1)+10*(j-1)+5*(q-1)+1}=[48*(i-1)+12*(j-1)+6*(q-1)+1,48*(i-1)+12*(j-1)+6*(q-1)+2,48*(i-1)+12*(j-1)+6*(q-1)+3];
            ori.newPanel{40*(i-1)+10*(j-1)+5*(q-1)+2}=[48*(i-1)+12*(j-1)+6*(q-1)+4,48*(i-1)+12*(j-1)+6*(q-1)+5,48*(i-1)+12*(j-1)+6*(q-1)+6];
            
            ori.newPanel{40*(i-1)+10*(j-1)+5*(q-1)+3}=[48*(i-1)+12*(j-1)+6*(q-1)+1,...
                48*(i-1)+12*(j-1)+6*(q-1)+4,...
                48*(i-1)+12*(j-1)+6*(q-1)+5,...
                48*(i-1)+12*(j-1)+6*(q-1)+2];
            
            ori.newPanel{40*(i-1)+10*(j-1)+5*(q-1)+4}=[48*(i-1)+12*(j-1)+6*(q-1)+1,...
                48*(i-1)+12*(j-1)+6*(q-1)+4,...
                48*(i-1)+12*(j-1)+6*(q-1)+6,...
                48*(i-1)+12*(j-1)+6*(q-1)+3];
            
            ori.newPanel{40*(i-1)+10*(j-1)+5*(q-1)+5}=[48*(i-1)+12*(j-1)+6*(q-1)+2,...
                48*(i-1)+12*(j-1)+6*(q-1)+5,...
                48*(i-1)+12*(j-1)+6*(q-1)+6,...
                48*(i-1)+12*(j-1)+6*(q-1)+3];
       
              
        end
    end
end

inPlateBarNum=96*N;

%% Here we start adding the connecting structures
for i=1:N
    
    for j=1:4
    
        if (mod(i,2)==1 && mod(j,2)==1) || (mod(i,2)==0 && mod(j,2)==0)
            
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+4, ...
                                      48*(i-1)+12*(j-1)+10,...
                                      48*(i-1)+12*(j-1)+6,...
                                      48*(i-1)+12*(j-1)+12,...
                                      48*(i-1)+12*(j-1)+1,...
                                      48*(i-1)+12*(j-1)+7,...
                                      48*(i-1)+12*(j-1)+3,...
                                      48*(i-1)+12*(j-1)+9,...
                                      connectorBarArea,connectorBarArea/5,L,gap,t,hingeStiff)
     
            if j~=4
                ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+11, ...
                          48*(i-1)+12*(j-1)+16,...
                          48*(i-1)+12*(j-1)+12,...
                          48*(i-1)+12*(j-1)+17,...
                          48*(i-1)+12*(j-1)+8,...
                          48*(i-1)+12*(j-1)+13,...
                          48*(i-1)+12*(j-1)+9,...
                          48*(i-1)+12*(j-1)+14,...
                          connectorBarArea,connectorBarArea/5,L,gap,t,hingeStiffYoshimura)
                
            end        
        else
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+6, ...
                          48*(i-1)+12*(j-1)+11,...
                          48*(i-1)+12*(j-1)+5,...
                          48*(i-1)+12*(j-1)+10,...
                          48*(i-1)+12*(j-1)+3,...
                          48*(i-1)+12*(j-1)+8,...
                          48*(i-1)+12*(j-1)+2,...
                          48*(i-1)+12*(j-1)+7,...
                          connectorBarArea,connectorBarArea/5,L,gap,t,hingeStiff)

            if j~=4
                
                ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+11, ...
                          48*(i-1)+12*(j-1)+16,...
                          48*(i-1)+12*(j-1)+12,...
                          48*(i-1)+12*(j-1)+17,...
                          48*(i-1)+12*(j-1)+8,...
                          48*(i-1)+12*(j-1)+13,...
                          48*(i-1)+12*(j-1)+9,...
                          48*(i-1)+12*(j-1)+14,...
                          connectorBarArea,connectorBarArea/5,L,gap,t,hingeStiffYoshimura)

            end
        end
    end    
end

%% Here we start adding the bottom connector
btmConnectorNum=inPlateBarNum+49*N;
for i=1:N-1
    for j=1:4
        if (mod(i,2)==1 && mod(j,2)==1) || (mod(i,2)==0 && mod(j,2)==0)
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+3, ...
              48*(i-1)+12*(j-1)+51,...
              48*(i-1)+12*(j-1)+2,...
              48*(i-1)+12*(j-1)+49,...
              48*(i-1)+12*(j-1)+6,...
              48*(i-1)+12*(j-1)+54,...
              48*(i-1)+12*(j-1)+5,...
              48*(i-1)+12*(j-1)+52,...
              connectorBarArea,connectorBarArea/5,L,gap,t,hingeStiff)
            
        else
            ori.AddHingeForThickPanel(48*(i-1)+12*(j-1)+7, ...
              48*(i-1)+12*(j-1)+55,...
              48*(i-1)+12*(j-1)+9,...
              48*(i-1)+12*(j-1)+56,...
              48*(i-1)+12*(j-1)+10,...
              48*(i-1)+12*(j-1)+58,...
              48*(i-1)+12*(j-1)+12,...
              48*(i-1)+12*(j-1)+59,...
              connectorBarArea,connectorBarArea/5,L,gap,t,hingeStiff)

        end
    end
end

    
%% Set up terms for initialization
newNodeNum = size(ori.newNode);
newNodeNum = newNodeNum(1);
ori.currentAppliedForce = zeros(newNodeNum,3);   
ori.currentU = zeros(newNodeNum,3);
ori.currentSprZeroStrain = pi*logical(ori.sprK);
end

