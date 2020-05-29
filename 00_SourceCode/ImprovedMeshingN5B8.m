%% Add additional structure for crease
% This code is used to generate an improved meshing for the bar and hinge
% model. Each hinge is further discretized into 6 more pannels and thus
% provide additional degree of freedom. The Panel is modeled with the
% improved N5B8 model.

% The content in this code is based on the bar and hinge model with 
% compliant creases developed by Y. Zhu and E. T. Filipov.
% [1] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 


%% Coding Info:
% Built: 2018 07 15 [Yi ZHU]
% Updated: 2019 01 22  [Yi ZHU]
%   This update changes the way additional crease structure is generated
%   This update also add a new BarType BarConnect and BarArea output
%   This update also add a new CreaseIJKL output for assembly of rotational
%   creases stiffness matrix
% Rewrite: 2018 07 30 [Yi ZHU] {For Nonlinear Solver}
% Rewrite: 2018 08 09 [Yi ZHU] {For N5B8 upgrade}
% Rewrite: 2019 01 22 [Yi ZHU] {Cleaning of code & better polygon support}

function[Node,Panel,BarType,BarConnect,BarArea,BarLength,...
    SprIJKL,SprTargetZeroStrain,SprK,Type1BarNum,OldCrease,...
    PanelInerBarStart,CenterNodeStart,NewFoldingSequence,...
    OldNode,PanelNum,oldPanel]...
    =ImprovedMeshingN5B8(Node0,Panel0,RotationZeroStrain,...
    FoldingSequence,ModelConstant)

    CreaseW=ModelConstant{1};
    PanelE=ModelConstant{2};
    CreaseE=ModelConstant{3};
    PanelThick=ModelConstant{4};
    CreaseThick=ModelConstant{5};
    PanelPoisson=ModelConstant{6};
    CreasePoisson=ModelConstant{7};
    Flag2D3D=ModelConstant{8};
    DiagonalRate=ModelConstant{9};
    LockingOpen=ModelConstant{10};
    CompliantCreaseOpen=ModelConstant{17};
    
    %% if Using Compliant Craese model with finite width
    if CompliantCreaseOpen ==1

        % Identify the creases of the old pattern for further processing
        [OldCreaseNum, OldCreaseConnect, OldCreaseType]=IdentifyCrease(Node0,Panel0);
        A=size(Node0);
        B=size(Panel0);
        PanelNum=B(2);
        NodeNum=A(1);
        Node=Node0;
        Panel=Panel0;

        count=1; 
        % Numbering system used in this code

        ZeroStrainPositionFactor=0.5;
        % This factor determines how the zero strain position is distributed
        % betweeen the three horizontal springs. 

        % These following two matrixes are introduced as a reference for generating the
        % additional structure of the system. 

        OldCrease=zeros(10,1);
        % This matrix stores coresponding old crease number of the new crease after
        % the panel is shrunk in its size

        OldNode=zeros(10,1);
        % This vector stores the coresponding old node number of the new node
        % after the panel is shrunk in its size    
        % It gives -1 if nodes are the center of panel or 0 if is in the middle
        % of the crease

        BarInfo=zeros(10,5);
        % This matrix stores how the creases are conected within the panel.
        % [1][2] are the two node of the crease, 
        % [3] is the one next to [1] while
        % [4] is the one next to [2];
        % [5] is the number of Panel;

        BarType=zeros(10,1);
        % This matrix stores the type of the bar after improving the meshing
        % [1] Panel Bars
        % [2] Vertical Crease Bars that connect different panels
        % [3] Diagonal Crease Bars
        % [4] Horizontal Crease Bars
        % [5] Panel Diagonal Bars
        
        BarConnect=zeros(10,2);
        % This matrix stoes how new bars are connected

        BarArea=zeros(10,1);
        % This matrix stores the Areas of bars
        BarLength=zeros(10,1);
        % This matrix stores the original length of bars

        SprIJKL=zeros(10,4);
        % This matrix stores the node information for rotational hinges
        SprTargetZeroStrain=zeros(10,1);
        % This matrix stores the zeros strain position of each creases
        SprK=zeros(10,1);
        % This matrix stores the rotation stiffness of each creases
        Type1BarNum=0;
        % number of type 1 bars
        
        NodalMass=zeros(10,1);
        % This Vector stores the mass assigned to each node
        
        % Note: I am abusing the MATLAB's ability to automaticly change the 
        % matrix size here. For huge origami sturctures, this part needs to
        % be reworked otherwise it could be extremely unefficient. But I do
        % not think that anybody will do the large scale computation with
        % this kind of matlab code anyway so I am happy to abuse the
        % ability. 
        
        NewFoldingSequence=zeros(10,1);
        % This vector stores which creases are folded first.

        %% shrink the panel size based on the prescribed creaseW
        for i=1:PanelNum
            B=Panel0(i);
            C=cell2mat(B);
            D=size(C);
            N=D(2);

            % Identify if the node need to be Converted
            % 0 do not need to be changed Pure Boundary Node
            % 1 need to be changed Partial Boundary Node
            % 2 need to be changed Internal Node
            TempNodeType=zeros(N,1);
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N));
            TempCrease(1,:)=[minNode maxNode];
            TempCreaseType=zeros(N,1);
            TempCreaseNum=zeros(N,1);

            for j=1:N-1
                minNode=min(C(j),C(j+1));
                maxNode=max(C(j),C(j+1));
                TempCrease(j+1,:)=[minNode maxNode];        
            end

            % Search for the Crease Type of each Crease
            for j=1:N
                for k=1:OldCreaseNum
                    if eq(TempCrease(j,:),OldCreaseConnect(k,:))
                        TempCreaseType(j)=OldCreaseType(k);
                        TempCreaseNum(j)=k;
                        break
                    end        
                end      
            end

            % Define the Node Type of each Node
            if TempCreaseType(1)>1
               TempNodeType(N)=TempNodeType(N)+1; 
            end    
            if TempCreaseType(N)>1
               TempNodeType(N)=TempNodeType(N)+1; 
            end

            for j=1:N-1
                if TempCreaseType(j)>1
                    TempNodeType(j)=TempNodeType(j)+1; 
                end
                if TempCreaseType(j+1)>1
                    TempNodeType(j)=TempNodeType(j)+1; 
                end   
            end

            %% Shrink the size of Panel in space
            %Calculate the node coordinate of panel  
            % N is the number of node inside one panel
            for j=1:N
                if TempNodeType(j)==0
                    Node(count,:)=Node0(C(j),:);  
                    OldNode(count)=C(j);                    
                    count=count+1;
                %Partial Boundary Node
                elseif TempNodeType(j)==1
                    V1=zeros(1,3);
                    V2=zeros(1,3);
                    Vector=zeros(1,3);
                    if j==1
                        V1=Node0(C(N),:)-Node0(C(1),:);
                        V2=Node0(C(2),:)-Node0(C(1),:);
                        V1=V1/norm(V1);
                        V2=V2/norm(V2);
                        if TempCreaseType(1)>1
                            Vector=V2;
                        elseif TempCreaseType(2)>1
                            Vector=V1;
                        end
                    elseif j==N
                        V1=Node0(C(1),:)-Node0(C(N),:);
                        V2=Node0(C(N-1),:)-Node0(C(N),:);
                        V1=V1/norm(V1);
                        V2=V2/norm(V2);    
                        if TempCreaseType(1)>1
                            Vector=V2;
                        elseif TempCreaseType(N)>1
                            Vector=V1;
                        end
                    else
                        V1=Node0(C(j+1),:)-Node0(C(j),:);
                        V2=Node0(C(j-1),:)-Node0(C(j),:);
                        V1=V1/norm(V1);
                        V2=V2/norm(V2);
                        if TempCreaseType(j)>1
                            Vector=V1;
                        elseif TempCreaseType(j+1)>1
                            Vector=V2;
                        end
                    end
                    if norm(cross(V1,V2))~=0
                        theta=acos(dot(V1,V2)/norm(V1)/norm(V2));
                        Length=abs(CreaseW/2/sin(pi-theta));
                        Vector=Length*Vector;
                        Node(count,:)=Node0(C(j),:)+Vector;
                        OldNode(count)=C(j);
                        count=count+1;
                    elseif norm(cross(V1,V2))==0
                        Node(count,:)=Node0(C(j),:);
                        OldNode(count)=C(j);
                        count=count+1;
                    end
                % Inner Node
                elseif TempNodeType(j)==2
                    V1=zeros(1,3);
                    V2=zeros(1,3);
                    if j==1
                        V1=Node0(C(N),:)-Node0(C(1),:);
                        V2=Node0(C(2),:)-Node0(C(1),:);
                        V1=V1/norm(V1);
                        V2=V2/norm(V2);
                    elseif j==N
                        V1=Node0(C(1),:)-Node0(C(N),:);
                        V2=Node0(C(N-1),:)-Node0(C(N),:);
                        V1=V1/norm(V1);
                        V2=V2/norm(V2);                
                    else
                        V1=Node0(C(j+1),:)-Node0(C(j),:);
                        V2=Node0(C(j-1),:)-Node0(C(j),:);
                        V1=V1/norm(V1);
                        V2=V2/norm(V2);
                    end
                    if norm(cross(V1,V2))~=0
                        theta=acos(dot(V1,V2)/norm(V1)/norm(V2));
                        Vector=V1+V2;
                        Vector=Vector/norm(Vector);
                        Vector=abs(CreaseW/2/sin(theta/2))*Vector;       
                        Node(count,:)=Node0(C(j),:)+Vector;
                        OldNode(count)=C(j);
                        count=count+1;
                    elseif norm(cross(V1,V2))==0
                        Node(count,:)=Node0(C(j),:)+Vector;
                        OldNode(count)=C(j);
                        count=count+1;
                    end
                end
            end

            % Formulate the connection information vector 
            % Formulate the bar length vector
            PanelOrder=zeros(1,N);    
            for j=1:N
               PanelOrder(j)=count-N-1+j;
               OldCrease(count-N-1+j)=TempCreaseNum(j);
               if j==1
                   BarInfo(count-N,:)=[count-N,count-1,count-N+1,count-2,i];
                   BarConnect(count-N,:)=[count-N,count-1];
                   BarLength(count-N)=norm(Node0(C(1),:)-Node0(C(N),:));          
               elseif j==2
                   BarInfo(count-N-1+j,:)=[count-N,count-N+1,count-1,count-N+2,i];
                   BarConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
                   BarLength(count-N-1+j)=norm(Node0(C(j-1),:)-Node0(C(j),:));
               elseif j==N
                   BarInfo(count-N-1+j,:)=[count-N-2+j,count-N-1+j,count-N-3+j,count-N,i];  
                   BarConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
                   BarLength(count-N-1+j)=norm(Node0(C(j-1),:)-Node0(C(j),:));
               else
                   BarInfo(count-N-1+j,:)=[count-N-2+j,count-N-1+j,count-N-3+j,count-N+j,i];
                   BarConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
                   BarLength(count-N-1+j)=norm(Node0(C(j-1),:)-Node0(C(j),:));
               end
            end

            % Assign Areas to bars of Panels (Peripheral)
            if N==2
            else        
                tempNode=zeros(1,3);
                for j=1:N
                    tempNode=tempNode+Node0(C(j),:);
                end
                centerNodetemp=tempNode/N;  
                % solve area Lsum
                area=0.5*norm(cross(Node0(C(1),:)-centerNodetemp,Node0(C(N),:)-centerNodetemp));
                Lsum=norm(centerNodetemp-Node0(C(1),:))+norm(Node0(C(N),:)-Node0(C(1),:));
                for j=2:N
                   area=area+0.5*norm(cross(Node0(C(j-1),:)-centerNodetemp,Node0(C(j),:)-centerNodetemp));
                   Lsum=Lsum+norm(centerNodetemp-Node0(C(j),:))+norm(Node0(C(j-1),:)-Node0(C(j),:));
                end

                for j=1:N
                    BarArea(count-N-1+j)=CalculateApanel(PanelThick(i),area,PanelPoisson,Lsum);
                end  
            end  
            Panel{i}=PanelOrder;    


        end
        BarType=ones(count-1,1);
        Type1BarNum=count-1;

        %% Generate the additional structure
        A=size(OldCrease);
        NewCreaseNum=A(1);
        panelCount=PanelNum+1;
        BarNum=count;

        for i=1:OldCreaseNum
            RotAtZeroStrain=(RotationZeroStrain(i)-pi)*ZeroStrainPositionFactor+pi;
            if OldCreaseType(i)==1
            else
                countTemp=0;
                edge1=0;
                edge2=0;
                for j=1:NewCreaseNum
                    if i==OldCrease(j)
                        if countTemp ==0
                            edge1=j;
                            countTemp=countTemp+1;
                        elseif countTemp==1
                            edge2=j;
                            countTemp=countTemp+1;
                        else
                            break
                        end                
                    end            
                end 
                NodeIndex11=BarInfo(edge1,1);
                NodeIndex12=BarInfo(edge1,2);
                NodeIndex21=BarInfo(edge2,1);
                NodeIndex22=BarInfo(edge2,2);  

                Node11=Node(BarInfo(edge1,1),:);
                Node12=Node(BarInfo(edge1,2),:);
                Node21=Node(BarInfo(edge2,1),:);
                Node22=Node(BarInfo(edge2,2),:);

                Node11Ref1=Node(BarInfo(edge1,3),:)-Node(BarInfo(edge1,1),:);
                Node11Ref2=Node(BarInfo(edge1,2),:)-Node(BarInfo(edge1,1),:);

                Node21Ref1=Node(BarInfo(edge1,4),:)-Node(BarInfo(edge1,2),:);
                Node21Ref2=Node(BarInfo(edge1,1),:)-Node(BarInfo(edge1,2),:);

                theta1=acos(dot(Node11Ref1,Node11Ref2)/norm(Node11Ref1)/norm(Node11Ref2));
                RefVector1=Node11Ref1/norm(Node11Ref1);
                RefVector1=RefVector1*CreaseW/2/sin(theta1);

                theta2=acos(dot(Node21Ref1,Node21Ref2)/norm(Node21Ref1)/norm(Node21Ref2));
                RefVector2=Node21Ref1/norm(Node21Ref1);
                RefVector2=RefVector2*CreaseW/2/sin(theta2);

                CreaseL=norm(Node11-Node12);

                % Check the right order of numbering if not correct, switch the
                % numbering oder automatically
                if norm(Node11-Node21)<norm(Node11-Node22)
                else
                    tempNode=Node22;
                    tempIndex=NodeIndex22;
                    Node22=Node21;
                    Node21=tempNode;
                    NodeIndex22=NodeIndex21;
                    NodeIndex21=tempIndex;
                end

                % Generate Three new node and their location for the analysis
                if Flag2D3D==3
                    Node(count,:)=Node11-RefVector1;
                    OldNode(count)=0;
                    count=count+1;
                    Node(count,:)=Node12-RefVector2;
                    OldNode(count)=0;
                    count=count+1;
                    Node(count,:)=0.5*(Node11+Node12)-0.5*(RefVector1+RefVector2);
                    OldNode(count)=0;
                    count=count+1;
                elseif Flag2D3D==2
                    Node(count,:)=(Node11+Node21)/2;
                    OldNode(count)=0;
                    count=count+1;
                    Node(count,:)=(Node12+Node22)/2;
                    OldNode(count)=0;
                    count=count+1;
                    Node(count,:)=0.25*(Node12+Node22+Node11+Node21);
                    OldNode(count)=0;
                    count=count+1;
                end

                % Generate New Panels added to the structure
                Panel{panelCount}=[count-1 count-3 NodeIndex11];
                oldPanel(panelCount)=0;
                panelCount=panelCount+1;
                Panel{panelCount}=[NodeIndex21 count-3 count-1];
                oldPanel(panelCount)=0;
                panelCount=panelCount+1;
                Panel{panelCount}=[NodeIndex21 count-1 NodeIndex22];
                oldPanel(panelCount)=0;
                panelCount=panelCount+1;
                Panel{panelCount}=[NodeIndex22 count-2 count-1];
                oldPanel(panelCount)=0;
                panelCount=panelCount+1;
                Panel{panelCount}=[NodeIndex12 count-1 count-2];
                oldPanel(panelCount)=0;
                panelCount=panelCount+1;
                Panel{panelCount}=[NodeIndex12 count-1 NodeIndex11];
                oldPanel(panelCount)=0;
                panelCount=panelCount+1;

                BarType(BarNum)=2;
                BarConnect(BarNum,:)=[NodeIndex11,count-3];
                BarArea(BarNum)=CalculateA1(CreaseL,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex11,:)-Node(count-3,:));
                SprIJKL(BarNum,:)=[0,0,0,0]; 
                BarNum=BarNum+1;

                BarType(BarNum)=3;
                BarConnect(BarNum,:)=[NodeIndex11,count-1];
                BarArea(BarNum)=CalculateA2(CreaseL,CreaseW,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex11,:)-Node(count-1,:));
                SprIJKL(BarNum,:)=[NodeIndex12,count-1,NodeIndex11,count-3];
                SprK(BarNum)=BarLength(BarNum)*CalculateKspr2(CreaseE,CreaseThick,CreaseW,DiagonalRate);
                NewFoldingSequence(BarNum)=FoldingSequence(i);
                SprTargetZeroStrain(BarNum)=pi;
                OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                BarType(BarNum)=3;
                BarConnect(BarNum,:)=[NodeIndex12,count-1];
                BarArea(BarNum)=CalculateA2(CreaseL,CreaseW,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex12,:)-Node(count-1,:));
                SprIJKL(BarNum,:)=[count-2,count-1,NodeIndex12,NodeIndex11];
                SprK(BarNum)=BarLength(BarNum)*CalculateKspr2(CreaseE,CreaseThick,CreaseW,DiagonalRate);
                NewFoldingSequence(BarNum)=FoldingSequence(i);
                SprTargetZeroStrain(BarNum)=pi;
                OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                BarType(BarNum)=2;
                BarConnect(BarNum,:)=[NodeIndex12,count-2];
                BarArea(BarNum)=CalculateA1(CreaseL,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex12,:)-Node(count-2,:));
                SprIJKL(BarNum,:)=[0,0,0,0];            
                BarNum=BarNum+1;

                BarType(BarNum)=4;
                BarConnect(BarNum,:)=[count-3,count-1];
                BarArea(BarNum)=CalculateA3(CreaseL,CreaseW,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(count-1,:)-Node(count-3,:));
                SprIJKL(BarNum,:)=[NodeIndex11,count-1,count-3,NodeIndex21];  
                SprK(BarNum)=BarLength(BarNum)*CalculateKspr1(CreaseE,CreaseThick,CreaseW);
                NewFoldingSequence(BarNum)=FoldingSequence(i);
                SprTargetZeroStrain(BarNum)=RotAtZeroStrain;
                OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                BarType(BarNum)=4;
                BarConnect(BarNum,:)=[count-2,count-1];
                BarArea(BarNum)=CalculateA3(CreaseL,CreaseW,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(count-2,:)-Node(count-1,:));
                SprIJKL(BarNum,:)=[NodeIndex12,count-2,count-1,NodeIndex22]; 
                SprK(BarNum)=BarLength(BarNum)*CalculateKspr1(CreaseE,CreaseThick,CreaseW);
                NewFoldingSequence(BarNum)=FoldingSequence(i);
                SprTargetZeroStrain(BarNum)=RotAtZeroStrain;
                OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                BarType(BarNum)=2;
                BarConnect(BarNum,:)=[NodeIndex21,count-3];
                BarArea(BarNum)=CalculateA1(CreaseL,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex21,:)-Node(count-3,:));
                SprIJKL(BarNum,:)=[0 0 0 0];
                BarNum=BarNum+1;

                BarType(BarNum)=3;
                BarConnect(BarNum,:)=[NodeIndex21,count-1];
                BarArea(BarNum)=CalculateA2(CreaseL,CreaseW,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex21,:)-Node(count-1,:));
                SprIJKL(BarNum,:)=[count-3,count-1,NodeIndex21,NodeIndex22];
                SprK(BarNum)=BarLength(BarNum)*CalculateKspr2(CreaseE,CreaseThick,CreaseW,DiagonalRate);
                NewFoldingSequence(BarNum)=FoldingSequence(i);
                SprTargetZeroStrain(BarNum)=pi;
                OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                BarType(BarNum)=3;
                BarConnect(BarNum,:)=[NodeIndex22,count-1];
                BarArea(BarNum)=CalculateA2(CreaseL,CreaseW,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex22,:)-Node(count-1,:));
                SprIJKL(BarNum,:)=[count-2,NodeIndex22,count-1,NodeIndex21];
                SprK(BarNum)=BarLength(BarNum)*CalculateKspr2(CreaseE,CreaseThick,CreaseW,DiagonalRate);
                NewFoldingSequence(BarNum)=FoldingSequence(i);
                SprTargetZeroStrain(BarNum)=pi;
                OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                BarType(BarNum)=2;
                BarConnect(BarNum,:)=[NodeIndex22,count-2];
                BarArea(BarNum)=CalculateA1(CreaseL,CreaseThick,CreasePoisson);
                BarLength(BarNum)=norm(Node(NodeIndex22,:)-Node(count-2,:));
                SprIJKL(BarNum,:)=[0 0 0 0];
                BarNum=BarNum+1;  
            end 
        end

        %% Add bars for the bending of panels
        CenterNodeStart=count-1;
        PanelInerBarStart=BarNum-1;
        for i=1:PanelNum
            B=Panel(i);
            C=cell2mat(B);
            D=size(C);
            N=D(2);

            if N==2
            else        
                tempNode=zeros(1,3);
                for j=1:N
                    tempNode=tempNode+Node(C(j),:);
                end
                Node(count,:)=tempNode/N;   
                OldNode(count)=-1;
                count=count+1;

                % solve area Lsum
                area=0.5*norm(cross(Node(C(1),:)-Node(count-1,:),Node(C(N),:)-Node(count-1,:)));
                Lsum=norm(Node(count-1,:)-Node(C(1),:))+norm(Node(C(N),:)-Node(C(1),:));
                for j=2:N
                   area=area+0.5*norm(cross(Node(C(j-1),:)-Node(count-1,:),Node(C(j),:)-Node(count-1,:)));
                   Lsum=Lsum+norm(Node(count-1,:)-Node(C(j),:))+norm(Node(C(j-1),:)-Node(C(j),:));
                end

                for j=1:N
                    BarType(BarNum)=5;
                    BarConnect(BarNum,:)=[C(j),count-1];
                    BarLength(BarNum)=norm(Node(C(j),:)-Node(count-1,:)); 
                    BarArea(BarNum)=CalculateApanel(PanelThick(i),area,PanelPoisson,Lsum); 
                    SprK(BarNum)=BarLength(BarNum)*CalculateKsprpanel(CreaseW,PanelThick(i),PanelE);

                    SprTargetZeroStrain(BarNum)=pi;
                    if j==1
                        SprIJKL(BarNum,:)=[C(N),C(1),count-1,C(2)];
                        Panel{panelCount}=[count-1,C(1),C(2)];
                        oldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    elseif j==N
                        SprIJKL(BarNum,:)=[C(N-1),C(N),count-1,C(1)];
                        Panel{panelCount}=[count-1,C(1),C(N)];
                        oldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    else
                        SprIJKL(BarNum,:)=[C(j-1),C(j),count-1,C(j+1)];
                        Panel{panelCount}=[count-1,C(j),C(j+1)];
                        oldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    end            
                    BarNum=BarNum+1;        
                end  
            end   
        end

        %% Calculate IJKL for panel bars
        count=1;
        for i=1:OldCreaseNum
            RotAtZeroStrain=(RotationZeroStrain(i)-pi)*(1-ZeroStrainPositionFactor)/2+pi;
            %RotAtZeroStrain=pi;
            if OldCreaseType(i)==1
            else
                countTemp=0;
                edge1=0;
                edge2=0;
                for j=1:NewCreaseNum
                    if i==OldCrease(j)
                        if countTemp ==0
                            edge1=j;
                            countTemp=countTemp+1;
                        elseif countTemp==1
                            edge2=j;
                            countTemp=countTemp+1;
                        else
                            break
                        end                
                    end            
                end 
                Node11=BarInfo(edge1,1);
                Node12=BarInfo(edge1,2);

                Node21=BarInfo(edge2,1);
                Node22=BarInfo(edge2,2);

                if norm(Node(Node11,:)-Node(Node21,:))<norm(Node(Node11,:)-Node(Node22,:))
                else
                    tempNode=Node22;
                    Node22=Node21;
                    Node21=tempNode;
                end

                Panel1=BarInfo(edge1,5);
                Panel2=BarInfo(edge2,5);
                CenterNode=Type1BarNum+3*(count-1)+3;

                % These are rotational springs at peripheral of panels, they use crease properties.
                SprIJKL(edge1,:)=[CenterNode,Node11,Node12,CenterNodeStart+Panel1];
                SprTargetZeroStrain(edge1)=RotAtZeroStrain;
                SprK(edge1)=BarLength(edge1)*CalculateKspr1(CreaseE,CreaseThick,CreaseW);
                NewFoldingSequence(edge1)=FoldingSequence(i);

                SprIJKL(edge2,:)=[CenterNode,Node22,Node21,CenterNodeStart+Panel2];
                SprTargetZeroStrain(edge2)=RotAtZeroStrain;
                SprK(edge2)=BarLength(edge2)*CalculateKspr1(CreaseE,CreaseThick,CreaseW);
                NewFoldingSequence(edge2)=FoldingSequence(i);

                count=count+1;
            end
        end
    %% If using the concentrated hinge model
    else
        A=size(Node0);
        B=size(Panel0);
        PanelNum=B(2);
        NodeNum=A(1);
        CenterNodeStart=NodeNum;
                
        Node=zeros(PanelNum+NodeNum,3);
        Node(1:NodeNum,1:3)=Node0(1:NodeNum,1:3);
        
        OldNode=zeros(PanelNum+NodeNum,1);
        for i=1:NodeNum
            OldNode(i)=i;
        end
        
        PanelCount=PanelNum+1;
        Panel=Panel0;
        
        BarCount=1;
        BarConnect=zeros(5,2);
        BarArea=zeros(5,1);
        BarLength=zeros(5,1);
        BarType=zeros(5,1);
        
        % This array stores how many time the bar appears, 
        % Value greater than 1 indicates The bar is shared by panels
        BarApperenceNum=ones(5,1);
        
        SprIJKL=zeros(5,4);
        SprK=zeros(5,1);
        SprTargetZeroStrain=zeros(5,1);
        NewFoldingSequence=zeros(5,1);        

        OldCrease=zeros(5,1);  
        
        % Fromulation of Panel and Node
        for i=1:PanelNum
            B=Panel0(i);
            C=cell2mat(B);
            N=size(C,2);     
            centerNode=zeros(1,3);           
            
            for j=1:N
                % Form the coordinates of center node
                centerNode=centerNode+Node(C(j),:);    
                % Form the panel information
                if j==1
                    Panel{PanelCount}=[NodeNum+i,C(1),C(N)];
                    oldPanel(panelCount)=i;
                else
                    Panel{PanelCount}=[NodeNum+i,C(j),C(j-1)];
                    oldPanel(panelCount)=i;
                end 
                PanelCount=PanelCount+1;
            end
            centerNode=centerNode/N;
            Node(i+NodeNum,:)=centerNode;
        end
        
        % Calculate Apanel
        Apanel=zeros(PanelNum,1);
        for i=1:PanelNum
            B=Panel0(i);
            C=cell2mat(B);
            N=size(C,2); 
            area=0;
            Lsum=0;
            
            for j=1:N                                
                if j==1
                    area=area+0.5*norm(cross(-Node(CenterNodeStart+i,:)+Node(C(N),:),-Node(CenterNodeStart+i,:)+Node(C(1),:)));
                    Lsum=Lsum+norm(-Node(CenterNodeStart+i,:)+Node(C(1),:))+norm(-Node(C(1),:)+Node(C(N),:));
                else
                    area=area+0.5*norm(cross(-Node(CenterNodeStart+i,:)+Node(C(j-1),:),-Node(CenterNodeStart+i,:)+Node(C(j),:)));
                    Lsum=Lsum+norm(-Node(CenterNodeStart+i,:)+Node(C(1),:))+norm(-Node(C(1),:)+Node(C(N),:));
                end
            end
            Apanel(i)=CalculateApanel(PanelThick(i),area,PanelPoisson,Lsum);
        end
        
        % Calculate information for type 1 bars (at the edge)
        for i=1:PanelNum
            B=Panel0(i);
            C=cell2mat(B);
            N=size(C,2);     
            centerNode=zeros(1,3);
            
            for j=1:N            
                if j==1
                    NodeMin=min(C(1),C(N));
                    NodeMax=max(C(1),C(N));
                    check=0;
                    BarIndex=0;
                    for k=1:BarCount-1
                       if (NodeMin==BarConnect(k,1) && NodeMax==BarConnect(k,2))
                           check=1;
                           BarIndex=k;
                           BarApperenceNum(BarIndex)=BarApperenceNum(BarIndex)+1;
                       end
                    end
                    if check==0
                       BarConnect(BarCount,:)=[NodeMin,NodeMax];
                       BarLength(BarCount,1)=norm(Node(NodeMin,:)-Node(NodeMax,:));
                       BarType(BarCount)=1;
                       BarArea(BarCount)=Apanel(i);
                       BarCount=BarCount+1;
                       BarApperenceNum(BarCount)=1;
                    else
                       BarArea(BarIndex)=BarArea(BarIndex)+Apanel(i);
                    end
                else
                    NodeMin=min(C(j),C(j-1));
                    NodeMax=max(C(j),C(j-1));
                    check=0;
                    BarIndex=0;
                    for k=1:BarCount-1
                       if (NodeMin==BarConnect(k,1) && NodeMax==BarConnect(k,2))
                           check=1;
                           BarIndex=k;
                           BarApperenceNum(BarIndex)=BarApperenceNum(BarIndex)+1;
                       end
                    end
                    if check==0
                       BarConnect(BarCount,:)=[NodeMin,NodeMax];
                       BarLength(BarCount,1)=norm(Node(NodeMin,:)-Node(NodeMax,:));
                       BarType(BarCount)=1;
                       BarArea(BarCount)=Apanel(i);
                       BarCount=BarCount+1;
                       BarApperenceNum(BarCount)=1;
                    else
                       BarArea(BarIndex)=BarArea(BarIndex)+Apanel(i);
                    end
                end
            end
        end
        
        PanelInerBarStart=BarCount-1;        
        Type1BarNum=BarCount-1;   

        % Calculate information for type 5 bars (bars inside the panels)
        for i=1:PanelNum
            B=Panel0(i);
            C=cell2mat(B);
            N=size(C,2);             
            Ksprpanel=CalculateKsprpanel(CreaseW,PanelThick(i),PanelE);
            
            for j=1:N  
                BarConnect(BarCount,:)=[CenterNodeStart+i,C(j)];
                BarLength(BarCount)=norm(Node(CenterNodeStart+i,:)-Node(C(j),:));                
                BarType(BarCount)=5;
                BarArea(BarCount)=Apanel(i);
                BarApperenceNum(BarCount)=1;                
                SprK(BarCount)=BarLength(BarCount)*Ksprpanel;
                SprTargetZeroStrain(BarCount)=pi;
                NewFoldingSequence(BarCount)=0;
                if j==1
                    SprIJKL(BarCount,:)=[C(2),C(1),NodeNum+i,C(N)];
                elseif j==N
                    SprIJKL(BarCount,:)=[C(1),C(N),NodeNum+i,C(N-1)];
                else
                    SprIJKL(BarCount,:)=[C(j+1),C(j),NodeNum+i,C(j-1)];                    
                end
                BarCount=BarCount+1;
            end
        end     
        
        % Formulate IJKL matrix for bars connecting different panels
        for i=1:BarCount-1
            if BarApperenceNum(i)>1
                nodeMin=BarConnect(i,1);
                nodeMax=BarConnect(i,2);
                temp=1;
                CenterIndex=zeros(2,1);
                
                for j=1:PanelNum
                     B=Panel0(j);
                     C=cell2mat(B);
                     N=size(C,2);
                     check=0;
                     for k=1:N
                         if nodeMin==C(k)
                             check=check+1;
                         end
                     end
                     for k=1:N
                         if nodeMax==C(k)
                             check=check+1;
                         end
                     end
                     
                     if check==2
                         CenterIndex(temp)=j;
                         temp=temp+1;                         
                     end
                end
                SprIJKL(i,:)=[NodeNum+CenterIndex(1),nodeMin,nodeMax,NodeNum+CenterIndex(2)];    
                SprK(i)=CalculateKhinge(CreaseW,CreaseThick,CreaseE)*BarLength(i);
                SprTargetZeroStrain(i)=RotationZeroStrain(i);
                NewFoldingSequence(i)=FoldingSequence(i);
            end
        end
    end
    Panel=Panel(PanelNum+1:end);
    oldPanel=oldPanel(PanelNum+1:end);
    
end


%% Functions used to define areas and rotational stiffness
function [A1]=CalculateA1(L,CreaseThick,CreasePoisson)
    A1=L*CreaseThick/2/(1-CreasePoisson^2);
end

function [A2]=CalculateA2(L,CreaseW,CreaseThick,CreasePoisson)
    A2=((L^2+CreaseW^2)^1.5)/4/(1+CreasePoisson)*CreaseThick/L/CreaseW;
end

function [A3]=CalculateA3(L,CreaseW,CreaseThick,CreasePoisson)
    A3=((L^2+CreaseW^2)^1.5)/4/(1+CreasePoisson)*CreaseThick/L/CreaseW;
end

function [Kspr1]=CalculateKspr1(Ecrease,CreaseThick,CreaseW)
    Kspr1=(Ecrease*CreaseThick^3)/4/CreaseW;
end

function [Kspr2]=CalculateKspr2(Ecrease,CreaseThick,CreaseW,DiagonalRate)
    Kspr2=DiagonalRate*(Ecrease*CreaseThick^3)/4/CreaseW;
end

function [Apanel]=CalculateApanel(PanelThick,area,PanelPoisson,Lsum)
    Apanel=2*PanelThick*area/(1-PanelPoisson)/Lsum;
end

function [Ksprpanel]=CalculateKsprpanel(CreaseW,PanelThick,PanelE)
    Ksprpanel=PanelE/12*(PanelThick^3)/CreaseW;
end

% This is for model with concentrated hinge
function [Khinge]=CalculateKhinge(CreaseW,CreaseThick,CreaseE)
    Khinge=CreaseE/12*(CreaseThick^3)/CreaseW;
end
