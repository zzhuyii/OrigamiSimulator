%% Gnerate the new geometry for origmai structure
% This function generate a new geometry for the origami structure 
% considering the existence of compliant creases;
%
% Input: 
%   [node0]: the original nodal coordinates;
%   {panel0}: the original panel informatoin;
%   oldCreaseNum: total number of the creases in the orignal system;
%   [oldCreaseConnect]: how the old creases are connected;
%   [oldCreaseType]: type of the old creaes: =1 for boundary crease; !=1 for middle crasese;
%
% Output:
%   [newNode]: the new nodal coordinates;
%   {newPanel}: the new panel information;
%   [barType]: the type of the new bars;
%       [1] Panel Bars
%       [2] Vertical Crease Bars that connect different panels
%       [3] Diagonal Crease Bars
%       [4] Horizontal Crease Bars
%       [5] Panel Diagonal Bars
%   [barConnect]: the two nodes that connects a bar;
%   [sprIJKL]: conentivity information of the rotaional springs;
%   type1BarNum: total number of panel bars;
%   panelInerBarStart: Number of the first type [5] bars;
%   centerNodeStart: Number of the first node at the center of paenl;
%   [newNode2OldNode]: mapping between the new node and the old node;
%   [newCrease2OldCrease]: mapping between the new crease and the old crease;
%   [newPanel2OldPanel]: mapping between the new panel and the old panel;
%   panelNum: total number of new panels;
%
%


function Mesh_CompliantCreaseGeometry(obj)

    %% if Using Compliant Craese model with finite width
    if obj.compliantCreaseOpen ==1
            
        % Identify the creases of the old pattern for further processing
        A=size(obj.node0);
        B=size(obj.panel0);
        obj.newPanelNum=B(2);
        nodeNum=A(1);
        obj.newNode=obj.node0;
        obj.newPanel=obj.panel0;

        count=1; 
        % Numbering system used in this code

        % These following two matrixes are introduced as a reference for 
        % generating the additional structure of the system. 

        obj.newCrease2OldCrease=zeros(10,1);
        % This matrix stores coresponding old crease number of the new 
        % crease after the panel is shrunk in its size

        obj.newNode2OldNode=zeros(10,1);
        % This vector stores the coresponding old node number of the new node
        % after the panel is shrunk in its size    
        % It gives -1 if nodes are the center of panel or 0 if is in the middle
        % of the crease

        barInfo=zeros(10,5);
        % This matrix stores how the creases are conected within the panel.
        % [1][2] are the two node of the crease, 
        % [3] is the one next to [1] while
        % [4] is the one next to [2];
        % [5] is the number of Panel;

        obj.barType=zeros(10,1);
        % This matrix stores the type of the bar after improving the meshing
        % [1] Panel Bars
        % [2] Vertical Crease Bars that connect different panels
        % [3] Diagonal Crease Bars
        % [4] Horizontal Crease Bars
        % [5] Panel Diagonal Bars
        
        obj.barConnect=zeros(10,2);
        % This matrix stoes how new bars are connected

        obj.sprIJKL=zeros(10,4);
        % This matrix stores the node information for rotational hinges
        
        obj.type1BarNum=0;
        % number of type 1 bars
        
%         obj.currentRotZeroStrain=pi*ones(obj,oldCreaseNum,1);
        

        %% shrink the panel size based on the prescribed creaseW
        for i=1:obj.newPanelNum
            
            % for the two loops, i is number of panel;
            % and j is the number of node inside a panel;
            % N is the total number of nodes within a panel;
            
            B=obj.panel0(i);
            C=cell2mat(B);
            D=size(C);
            N=D(2);

            % Identify if the node need to be Converted
            % 0 do not need to be changed Pure Boundary Node
            % 1 need to be changed Partial Boundary Node
            % 2 need to be changed Internal Node
            tempNodeType=zeros(N,1);
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N));
            tempCrease(1,:)=[minNode maxNode];
            tempCreaseType=zeros(N,1);
            tempCreaseNum=zeros(N,1);

            for j=1:N-1
                minNode=min(C(j),C(j+1));
                maxNode=max(C(j),C(j+1));
                tempCrease(j+1,:)=[minNode maxNode];        
            end

            % Search for the Crease Type of each Crease
            for j=1:N
                for k=1:obj.oldCreaseNum
                    if eq(tempCrease(j,:),obj.oldCreaseConnect(k,:))
                        tempCreaseType(j)=obj.oldCreaseType(k);
                        tempCreaseNum(j)=k;
                        break
                    end        
                end      
            end

            % Define the Node Type of each Node
            if tempCreaseType(1)>1
               tempNodeType(N)=tempNodeType(N)+1; 
            end    
            if tempCreaseType(N)>1
               tempNodeType(N)=tempNodeType(N)+1; 
            end

            for j=1:N-1
                if tempCreaseType(j)>1
                    tempNodeType(j)=tempNodeType(j)+1; 
                end
                if tempCreaseType(j+1)>1
                    tempNodeType(j)=tempNodeType(j)+1; 
                end   
            end

            %% Shrink the size of Panel in space
            % Calculate the node coordinate of panel  
            % N is the number of node inside one panel
            for j=1:N
                if tempNodeType(j)==0
                    obj.newNode(count,:)=obj.node0(C(j),:);  
                    obj.newNode2OldNode(count)=C(j);                    
                    count=count+1;
                % Partial Boundary Node
                elseif tempNodeType(j)==1
                    v1=zeros(1,3);
                    v2=zeros(1,3);
                    vector=zeros(1,3);
                    creaseW=0;
                    if j==1
                        v1=obj.node0(C(N),:)-obj.node0(C(1),:);
                        v2=obj.node0(C(2),:)-obj.node0(C(1),:);
                        
                        crease1Num=SearchCreaseNum(obj.oldCreaseConnect,C(1),C(N));
                        crease2Num=SearchCreaseNum(obj.oldCreaseConnect,C(2),C(1));
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        if tempCreaseType(1)>1
                            vector=v2;
                            creaseW=obj.creaseWidthVec(crease1Num);
                        elseif tempCreaseType(2)>1
                            vector=v1;
                            creaseW=obj.creaseWidthVec(crease2Num);
                        end
                    elseif j==N
                        v1=obj.node0(C(1),:)-obj.node0(C(N),:);
                        v2=obj.node0(C(N-1),:)-obj.node0(C(N),:);
                        
                        crease1Num=SearchCreaseNum(obj.oldCreaseConnect,C(1),C(N));
                        crease2Num=SearchCreaseNum(obj.oldCreaseConnect,C(N-1),C(N));
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);    
                        if tempCreaseType(1)>1
                            vector=v2;
                            creaseW=obj.creaseWidthVec(crease1Num);
                        elseif tempCreaseType(N)>1
                            vector=v1;
                            creaseW=obj.creaseWidthVec(crease2Num);
                        end
                    else
                        v1=obj.node0(C(j+1),:)-obj.node0(C(j),:);
                        v2=obj.node0(C(j-1),:)-obj.node0(C(j),:);
                        
                        crease1Num=SearchCreaseNum(obj.oldCreaseConnect,C(j+1),C(j));
                        crease2Num=SearchCreaseNum(obj.oldCreaseConnect,C(j-1),C(j));
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        if tempCreaseType(j)>1
                            vector=v1;
                            creaseW=obj.creaseWidthVec(crease2Num);
                        elseif tempCreaseType(j+1)>1
                            vector=v2;
                            creaseW=obj.creaseWidthVec(crease1Num);
                        end
                    end
                    if norm(cross(v1,v2))~=0
                        theta=acos(dot(v1,v2)/norm(v1)/norm(v2));
                        Length=abs(creaseW/2/sin(pi-theta));
                        vector=Length*vector;
                        obj.newNode(count,:)=obj.node0(C(j),:)+vector;
                        obj.newNode2OldNode(count)=C(j);
                        count=count+1;
                    elseif norm(cross(v1,v2))==0
                        obj.newNode(count,:)=obj.node0(C(j),:);
                        obj.newNode2OldNode(count)=C(j);
                        count=count+1;
                    end
                % Inner Node
                elseif tempNodeType(j)==2
                    v1=zeros(1,3);
                    v2=zeros(1,3);
                    if j==1
                        v1=obj.node0(C(N),:)-obj.node0(C(1),:);
                        v2=obj.node0(C(2),:)-obj.node0(C(1),:);
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        crease1Num=SearchCreaseNum(obj.oldCreaseConnect,C(N),C(1));
                        crease2Num=SearchCreaseNum(obj.oldCreaseConnect,C(2),C(1));
                    elseif j==N
                        v1=obj.node0(C(1),:)-obj.node0(C(N),:);
                        v2=obj.node0(C(N-1),:)-obj.node0(C(N),:);
                        v1=v1/norm(v1);
                        v2=v2/norm(v2); 
                        crease1Num=SearchCreaseNum(obj.oldCreaseConnect,C(1),C(N));
                        crease2Num=SearchCreaseNum(obj.oldCreaseConnect,C(N-1),C(N));
                    else
                        v1=obj.node0(C(j+1),:)-obj.node0(C(j),:);
                        v2=obj.node0(C(j-1),:)-obj.node0(C(j),:);
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        crease1Num=SearchCreaseNum(obj.oldCreaseConnect,C(j+1),C(j));
                        crease2Num=SearchCreaseNum(obj.oldCreaseConnect,C(j-1),C(j));
                    end
                    if norm(cross(v1,v2))~=0
                        crease1W=obj.creaseWidthVec(crease1Num);
                        crease2W=obj.creaseWidthVec(crease2Num);
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        
                        vector=crease1W/2*v2+crease2W/2*v1;
                        
                        obj.newNode(count,:)=obj.node0(C(j),:)+vector;
                        obj.newNode2OldNode(count)=C(j);
                        count=count+1;
                    elseif norm(cross(v1,v2))==0
                        obj.newNode(count,:)=obj.node0(C(j),:)+vector;
                        obj.newNode2OldNode(count)=C(j);
                        count=count+1;
                    end
                end
            end

            % Formulate the connection information vector 
            % Formulate the bar length vector
            PanelOrder=zeros(1,N);    
            for j=1:N
               PanelOrder(j)=count-N-1+j;
               obj.newCrease2OldCrease(count-N-1+j)=tempCreaseNum(j);
               if j==1
                   barInfo(count-N,:)=[count-N,count-1,count-N+1,count-2,i];
                   obj.barConnect(count-N,:)=[count-N,count-1];      
               elseif j==2
                   barInfo(count-N-1+j,:)=[count-N,count-N+1,count-1,count-N+2,i];
                   obj.barConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
               elseif j==N
                   barInfo(count-N-1+j,:)=[count-N-2+j,count-N-1+j,count-N-3+j,count-N,i];  
                   obj.barConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
               else
                   barInfo(count-N-1+j,:)=[count-N-2+j,count-N-1+j,count-N-3+j,count-N+j,i];
                   obj.barConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
               end
            end

            % Assign Areas to bars of Panels (Peripheral)
            if N==2
            else        
                tempNode=zeros(1,3);
                for j=1:N
                    tempNode=tempNode+obj.node0(C(j),:);
                end
                centerNodetemp=tempNode/N;  
            end  
            obj.newPanel{i}=PanelOrder;    


        end
        obj.barType=ones(count-1,1);
        obj.type1BarNum=count-1;

        %% Generate the additional structure
        A=size(obj.newCrease2OldCrease);
        NewCreaseNum=A(1);
        panelCount=obj.newPanelNum+1;
        BarNum=count;

        for i=1:obj.oldCreaseNum
            creaseW = obj.creaseWidthVec(i);
            if obj.oldCreaseType(i)==1
            else
                countTemp=0;
                edge1=0;
                edge2=0;
                for j=1:NewCreaseNum
                    if i==obj.newCrease2OldCrease(j)
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
                NodeIndex11=barInfo(edge1,1);
                NodeIndex12=barInfo(edge1,2);
                NodeIndex21=barInfo(edge2,1);
                NodeIndex22=barInfo(edge2,2);  

                Node11=obj.newNode(barInfo(edge1,1),:);
                Node12=obj.newNode(barInfo(edge1,2),:);
                Node21=obj.newNode(barInfo(edge2,1),:);
                Node22=obj.newNode(barInfo(edge2,2),:);

                Node11Ref1=obj.newNode(barInfo(edge1,3),:)-obj.newNode(barInfo(edge1,1),:);
                Node11Ref2=obj.newNode(barInfo(edge1,2),:)-obj.newNode(barInfo(edge1,1),:);

                Node21Ref1=obj.newNode(barInfo(edge1,4),:)-obj.newNode(barInfo(edge1,2),:);
                Node21Ref2=obj.newNode(barInfo(edge1,1),:)-obj.newNode(barInfo(edge1,2),:);

                theta1=acos(dot(Node11Ref1,Node11Ref2)/norm(Node11Ref1)/norm(Node11Ref2));
                RefVector1=Node11Ref1/norm(Node11Ref1);
                RefVector1=RefVector1*creaseW/2/sin(theta1);

                theta2=acos(dot(Node21Ref1,Node21Ref2)/norm(Node21Ref1)/norm(Node21Ref2));
                RefVector2=Node21Ref1/norm(Node21Ref1);
                RefVector2=RefVector2*creaseW/2/sin(theta2);

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
                
                % calculate the current RotZero
%                 paenlVec11 = obj.newNode(barInfo(edge1,3),:)-obj.newNode(barInfo(edge1,1),:);
%                 paenlVec12 = obj.newNode(barInfo(edge1,2),:)-obj.newNode(barInfo(edge1,1),:);                
%                 paenlVec21 = obj.newNode(barInfo(edge2,3),:)-obj.newNode(barInfo(edge2,1),:);
%                 paenlVec22 = obj.newNode(barInfo(edge2,2),:)-obj.newNode(barInfo(edge2,1),:);                
%                 normalVec1=cross(paenlVec11/norm(paenlVec11),paenlVec12/norm(paenlVec12));
%                 normalVec2=cross(paenlVec21/norm(paenlVec21),paenlVec22/norm(paenlVec22));
%                 obj.currentRotZeroStrain(i)=

                % Generate Three new node and their location for the analysis
                if obj.mesh2D3D==3
%                     p1=(Node11-RefVector1)...
%                         *((sin(theta1/2)-sin(theta1/4))/sin(theta1/2))...
%                         +((Node11+Node21)/2)...
%                         *(sin(theta1/4)/sin(theta1/2));
                    p1=(Node11-RefVector1)...
                        *1/2 ...
                        +((Node11+Node21)/2)...
                        *1/2;
                    obj.newNode2OldNode(count)=0;
                    obj.newNode(count,:)=p1;
                    count=count+1;
                    
%                     p2=(Node12-RefVector2)...
%                         *((sin(theta2/2)-sin(theta2/4))/sin(theta2/2))...
%                         +((Node12+Node22)/2)...
%                         *(sin(theta2/4)/sin(theta2/2));
                    p2=(Node12-RefVector2)...
                        *1/2 ...
                        +((Node12+Node22)/2)...
                        *1/2;
                    obj.newNode2OldNode(count)=0;
                    obj.newNode(count,:)=p2;
                    count=count+1;
                    
                    obj.newNode(count,:)=(p1+p2)/2;
                    obj.newNode2OldNode(count)=0;
                    count=count+1;
                    
                elseif obj.mesh2D3D==2
                    obj.newNode(count,:)=(Node11+Node21)/2;
                    obj.newNode2OldNode(count)=0;
                    count=count+1;
                    obj.newNode(count,:)=(Node12+Node22)/2;
                    obj.newNode2OldNode(count)=0;
                    count=count+1;
                    obj.newNode(count,:)=0.25*(Node12+Node22+Node11+Node21);
                    obj.newNode2OldNode(count)=0;
                    count=count+1;
                end

                % Generate New Panels added to the structure
                obj.newPanel{panelCount}=[count-1 count-3 NodeIndex11];
                obj.newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                obj.newPanel{panelCount}=[NodeIndex21 count-3 count-1];
                obj.newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                obj.newPanel{panelCount}=[NodeIndex21 count-1 NodeIndex22];
                obj.newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                obj.newPanel{panelCount}=[NodeIndex22 count-2 count-1];
                obj.newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                obj.newPanel{panelCount}=[NodeIndex12 count-1 count-2];
                obj.newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                obj.newPanel{panelCount}=[NodeIndex12 count-1 NodeIndex11];
                obj.newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;

                obj.barType(BarNum)=2;
                obj.barConnect(BarNum,:)=[NodeIndex11,count-3];
                obj.sprIJKL(BarNum,:)=[0,0,0,0]; 
                BarNum=BarNum+1;

                obj.barType(BarNum)=3;
                obj.barConnect(BarNum,:)=[NodeIndex11,count-1];                
                obj.sprIJKL(BarNum,:)=[NodeIndex12,count-1,NodeIndex11,count-3];                
                obj.newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                obj.barType(BarNum)=3;
                obj.barConnect(BarNum,:)=[NodeIndex12,count-1];
                obj.sprIJKL(BarNum,:)=[count-2,count-1,NodeIndex12,NodeIndex11];
                obj.newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                obj.barType(BarNum)=2;
                obj.barConnect(BarNum,:)=[NodeIndex12,count-2];
                obj.sprIJKL(BarNum,:)=[0,0,0,0];            
                BarNum=BarNum+1;

                obj.barType(BarNum)=4;
                obj.barConnect(BarNum,:)=[count-3,count-1];
                obj.sprIJKL(BarNum,:)=[NodeIndex11,count-1,count-3,NodeIndex21];  
                obj.newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                obj.barType(BarNum)=4;
                obj.barConnect(BarNum,:)=[count-2,count-1];
                obj.sprIJKL(BarNum,:)=[NodeIndex12,count-2,count-1,NodeIndex22]; 
                obj.newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                obj.barType(BarNum)=2;
                obj.barConnect(BarNum,:)=[NodeIndex21,count-3];
                obj.sprIJKL(BarNum,:)=[0 0 0 0];
                BarNum=BarNum+1;

                obj.barType(BarNum)=3;
                obj.barConnect(BarNum,:)=[NodeIndex21,count-1];
                obj.sprIJKL(BarNum,:)=[count-3,count-1,NodeIndex21,NodeIndex22];
                obj.newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                obj.barType(BarNum)=3;
                obj.barConnect(BarNum,:)=[NodeIndex22,count-1];
                obj.sprIJKL(BarNum,:)=[count-2,NodeIndex22,count-1,NodeIndex21];
                obj.newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                obj.barType(BarNum)=2;
                obj.barConnect(BarNum,:)=[NodeIndex22,count-2];
                obj.sprIJKL(BarNum,:)=[0 0 0 0];
                BarNum=BarNum+1;  
            end 
        end

        %% Add bars for the bending of panels
        obj.centerNodeStart=count-1;
        obj.panelInnerBarStart=BarNum-1;
        for i=1:obj.newPanelNum
            B=obj.newPanel(i);
            C=cell2mat(B);
            D=size(C);
            N=D(2);

            if N==2
            else        
                tempNode=zeros(1,3);
                for j=1:N
                    tempNode=tempNode+obj.newNode(C(j),:);
                end
                obj.newNode(count,:)=tempNode/N;   
                obj.newNode2OldNode(count)=-1;
                count=count+1;

                for j=1:N
                    obj.barType(BarNum)=5;
                    obj.barConnect(BarNum,:)=[C(j),count-1];
                    if j==1
                        obj.sprIJKL(BarNum,:)=[C(N),C(1),count-1,C(2)];
                        obj.newPanel{panelCount}=[count-1,C(1),C(2)];
                        obj.newPanel2OldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    elseif j==N
                        obj.sprIJKL(BarNum,:)=[C(N-1),C(N),count-1,C(1)];
                        obj.newPanel{panelCount}=[count-1,C(1),C(N)];
                        obj.newPanel2OldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    else
                        obj.sprIJKL(BarNum,:)=[C(j-1),C(j),count-1,C(j+1)];
                        obj.newPanel{panelCount}=[count-1,C(j),C(j+1)];
                        obj.newPanel2OldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    end            
                    BarNum=BarNum+1;        
                end  
            end   
        end

        %% Calculate IJKL for panel bars
        count=1;
        for i=1:obj.oldCreaseNum
            if obj.oldCreaseType(i)==1
            else
                countTemp=0;
                edge1=0;
                edge2=0;
                for j=1:NewCreaseNum
                    if i==obj.newCrease2OldCrease(j)
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
                Node11=barInfo(edge1,1);
                Node12=barInfo(edge1,2);

                Node21=barInfo(edge2,1);
                Node22=barInfo(edge2,2);

                if norm(obj.newNode(Node11,:)-obj.newNode(Node21,:))<norm(obj.newNode(Node11,:)-obj.newNode(Node22,:))
                else
                    tempNode=Node22;
                    Node22=Node21;
                    Node21=tempNode;
                end

                Panel1=barInfo(edge1,5);
                Panel2=barInfo(edge2,5);
                CenterNode=obj.type1BarNum+3*(count-1)+3;

                % These are rotational springs at peripheral of panels, they use crease properties.
                obj.sprIJKL(edge1,:)=[CenterNode,Node11,Node12,obj.centerNodeStart+Panel1];
                obj.sprIJKL(edge2,:)=[CenterNode,Node22,Node21,obj.centerNodeStart+Panel2];
                count=count+1;
            end
        end               

        
    %% If using the concentrated hinge model         
    else 

        A=size(obj.node0);
        B=size(obj.panel0);
        obj.newPanelNum=B(2);
        nodeNum=A(1);
        obj.centerNodeStart=nodeNum;
                
        obj.newNode=zeros(obj.newPanelNum+nodeNum,3);
        obj.newNode(1:nodeNum,1:3)=obj.node0(1:nodeNum,1:3);
        
        obj.newNode2OldNode=zeros(obj.newPanelNum+nodeNum,1);
        for i=1:nodeNum
            obj.newNode2OldNode(i)=i;
        end
        
        panelCount=obj.newPanelNum+1;
        obj.newPanel=obj.panel0;
        
        barCount=1;
        obj.barConnect=zeros(5,2);
        obj.barType=zeros(5,1);
        
        % This array stores how many time the bar appears, 
        % Value greater than 1 indicates The bar is shared by panels
        barApperenceNum=ones(5,1);
        
        obj.sprIJKL=zeros(5,4);
        obj.newCrease2OldCrease=zeros(5,1);  
        
        % Fromulation of Panel and Node
        for i=1:obj.newPanelNum
            B=obj.panel0(i);
            C=cell2mat(B);
            N=size(C,2);     
            centerNode=zeros(1,3);           
            
            for j=1:N
                % Form the coordinates of center node
                centerNode=centerNode+obj.newNode(C(j),:);    
                % Form the panel information
                if j==1
                    obj.newPanel{panelCount}=[nodeNum+i,C(1),C(N)];
                    obj.newPanel2OldPanel(panelCount)=i;
                else
                    obj.newPanel{panelCount}=[nodeNum+i,C(j),C(j-1)];
                    obj.newPanel2OldPanel(panelCount)=i;
                end 
                panelCount=panelCount+1;
            end
            centerNode=centerNode/N;
            obj.newNode(i+nodeNum,:)=centerNode;
        end
        
        % Calculate information for type 1 bars (at the edge)
        for i=1:obj.newPanelNum
            B=obj.panel0(i);
            C=cell2mat(B);
            N=size(C,2);     
            centerNode=zeros(1,3);
            
            for j=1:N            
                if j==1
                    NodeMin=min(C(1),C(N));
                    NodeMax=max(C(1),C(N));
                    check=0;
                    barIndex=0;
                    for k=1:barCount-1
                       if (NodeMin==obj.barConnect(k,1) && NodeMax==obj.barConnect(k,2))
                           check=1;
                           barIndex=k;
                           barApperenceNum(barIndex)=barApperenceNum(barIndex)+1;
                       end
                    end
                    if check==0
                       obj.barConnect(barCount,:)=[NodeMin,NodeMax];
                       obj.barType(barCount)=1;
                       barCount=barCount+1;
                       barApperenceNum(barCount)=1;
                    else
                    end
                else
                    NodeMin=min(C(j),C(j-1));
                    NodeMax=max(C(j),C(j-1));
                    check=0;
                    barIndex=0;
                    for k=1:barCount-1
                       if (NodeMin==obj.barConnect(k,1) && NodeMax==obj.barConnect(k,2))
                           check=1;
                           barIndex=k;
                           barApperenceNum(barIndex)=barApperenceNum(barIndex)+1;
                       end
                    end
                    if check==0
                       obj.barConnect(barCount,:)=[NodeMin,NodeMax];
                       obj.barType(barCount)=1;
                       barCount=barCount+1;
                       barApperenceNum(barCount)=1;
                    else
                    end
                end
            end
        end
        
        obj.panelInnerBarStart=barCount-1;        
        obj.type1BarNum=barCount-1;   

        % Calculate information for type 5 bars (bars inside the panels)
        for i=1:obj.newPanelNum
            B=obj.panel0(i);
            C=cell2mat(B);
            N=size(C,2);             
            
            for j=1:N  
                obj.barConnect(barCount,:)=[C(j),obj.centerNodeStart+i];               
                obj.barType(barCount)=5;
                barApperenceNum(barCount)=1;       
                if j==1
                    obj.sprIJKL(barCount,:)=[C(2),C(1),nodeNum+i,C(N)];
                elseif j==N
                    obj.sprIJKL(barCount,:)=[C(1),C(N),nodeNum+i,C(N-1)];
                else
                    obj.sprIJKL(barCount,:)=[C(j+1),C(j),nodeNum+i,C(j-1)];                    
                end
                barCount=barCount+1;
            end
        end     
        
        % Formulate IJKL matrix for bars connecting different panels
        for i=1:barCount-1
            if barApperenceNum(i)>1
                nodeMin=obj.barConnect(i,1);
                nodeMax=obj.barConnect(i,2);
                temp=1;
                CenterIndex=zeros(2,1);
                
                for j=1:obj.newPanelNum
                     B=obj.panel0(j);
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
                obj.sprIJKL(i,:)=[nodeNum+CenterIndex(1),nodeMin,nodeMax,nodeNum+CenterIndex(2)];    
            end
        end
    end
    obj.newPanel=obj.newPanel(obj.newPanelNum+1:end);
    obj.newPanel2OldPanel=obj.newPanel2OldPanel(obj.newPanelNum+1:end);
    
    A=size(obj.barConnect);
    barNum=A(1);
    
    %% check the IJKL matrix to make sure the right polarity
    for t=1:barNum
       if obj.sprIJKL(t,1)~=0
           ni=obj.newNode(obj.sprIJKL(t,1),:);
           nj=obj.newNode(obj.sprIJKL(t,2),:);
           nk=obj.newNode(obj.sprIJKL(t,3),:);
           
           rij=nj-ni;
           rkj=nk-nj;
           
           vec=cross(rij,rkj);
           if vec(3)>0
               tempIndex=obj.sprIJKL(t,1);
               obj.sprIJKL(t,1)=obj.sprIJKL(t,4);
               obj.sprIJKL(t,4)=tempIndex;
           end
       end
    end
    
    
    %% finalize result
    obj.newCrease2OldCrease(obj.panelInnerBarStart:barNum)=zeros(size(obj.panelInnerBarStart:barNum));
    
    obj.Mesh_NumberingForContact()
    obj.Mesh_BarLength()
end

function creaseNum=SearchCreaseNum(oldCreaseConnect,i,j)
    minIndex=min(i,j);
    maxIndex=max(i,j);
    A=size(oldCreaseConnect);
    N=A(1);
    creaseNum=0;
    for i=1:N
        if oldCreaseConnect(i,1)==minIndex && oldCreaseConnect(i,2)==maxIndex
            creaseNum=i;
        end
    end
end