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


function[newNode,newPanel,barType,barConnect,...
    sprIJKL,type1BarNum,panelInnerBarStart,centerNodeStart,...
    newNode2OldNode,newCrease2OldCrease,newPanel2OldPanel,newPanelNum]...
    =Mesh_CompliantCreaseGeometry(node0,panel0,...
    oldCreaseNum,oldCreaseConnect,oldCreaseType,...
    modelGeometryConstant)

    
    Flag2D3D=modelGeometryConstant{1};
    CompliantCreaseOpen=modelGeometryConstant{2};
    creaseWidthMat=modelGeometryConstant{3};
    
    
    %% if Using Compliant Craese model with finite width
    if CompliantCreaseOpen ==1
            
        % Identify the creases of the old pattern for further processing
        A=size(node0);
        B=size(panel0);
        newPanelNum=B(2);
        nodeNum=A(1);
        newNode=node0;
        newPanel=panel0;

        count=1; 
        % Numbering system used in this code

        % These following two matrixes are introduced as a reference for generating the
        % additional structure of the system. 

        newCrease2OldCrease=zeros(10,1);
        % This matrix stores coresponding old crease number of the new crease after
        % the panel is shrunk in its size

        newNode2OldNode=zeros(10,1);
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

        barType=zeros(10,1);
        % This matrix stores the type of the bar after improving the meshing
        % [1] Panel Bars
        % [2] Vertical Crease Bars that connect different panels
        % [3] Diagonal Crease Bars
        % [4] Horizontal Crease Bars
        % [5] Panel Diagonal Bars
        
        barConnect=zeros(10,2);
        % This matrix stoes how new bars are connected

        sprIJKL=zeros(10,4);
        % This matrix stores the node information for rotational hinges
        
        type1BarNum=0;
        % number of type 1 bars
        

        %% shrink the panel size based on the prescribed creaseW
        for i=1:newPanelNum
            
            % for the two loops, i is number of panel;
            % and j is the number of node inside a panel;
            % N is the total number of nodes within a panel;
            
            B=panel0(i);
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
                for k=1:oldCreaseNum
                    if eq(tempCrease(j,:),oldCreaseConnect(k,:))
                        tempCreaseType(j)=oldCreaseType(k);
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
                    newNode(count,:)=node0(C(j),:);  
                    newNode2OldNode(count)=C(j);                    
                    count=count+1;
                % Partial Boundary Node
                elseif tempNodeType(j)==1
                    v1=zeros(1,3);
                    v2=zeros(1,3);
                    vector=zeros(1,3);
                    creaseW=0;
                    if j==1
                        v1=node0(C(N),:)-node0(C(1),:);
                        v2=node0(C(2),:)-node0(C(1),:);
                        
                        crease1Num=SearchCreaseNum(oldCreaseConnect,C(1),C(N));
                        crease2Num=SearchCreaseNum(oldCreaseConnect,C(2),C(1));
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        if tempCreaseType(1)>1
                            vector=v2;
                            creaseW=creaseWidthMat(crease1Num);
                        elseif tempCreaseType(2)>1
                            vector=v1;
                            creaseW=creaseWidthMat(crease2Num);
                        end
                    elseif j==N
                        v1=node0(C(1),:)-node0(C(N),:);
                        v2=node0(C(N-1),:)-node0(C(N),:);
                        
                        crease1Num=SearchCreaseNum(oldCreaseConnect,C(1),C(N));
                        crease2Num=SearchCreaseNum(oldCreaseConnect,C(N-1),C(N));
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);    
                        if tempCreaseType(1)>1
                            vector=v2;
                            creaseW=creaseWidthMat(crease1Num);
                        elseif tempCreaseType(N)>1
                            vector=v1;
                            creaseW=creaseWidthMat(crease2Num);
                        end
                    else
                        v1=node0(C(j+1),:)-node0(C(j),:);
                        v2=node0(C(j-1),:)-node0(C(j),:);
                        
                        crease1Num=SearchCreaseNum(oldCreaseConnect,C(j+1),C(j));
                        crease2Num=SearchCreaseNum(oldCreaseConnect,C(j-1),C(j));
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        if tempCreaseType(j)>1
                            vector=v1;
                            creaseW=creaseWidthMat(crease2Num);
                        elseif tempCreaseType(j+1)>1
                            vector=v2;
                            creaseW=creaseWidthMat(crease1Num);
                        end
                    end
                    if norm(cross(v1,v2))~=0
                        theta=acos(dot(v1,v2)/norm(v1)/norm(v2));
                        Length=abs(creaseW/2/sin(pi-theta));
                        vector=Length*vector;
                        newNode(count,:)=node0(C(j),:)+vector;
                        newNode2OldNode(count)=C(j);
                        count=count+1;
                    elseif norm(cross(v1,v2))==0
                        newNode(count,:)=node0(C(j),:);
                        newNode2OldNode(count)=C(j);
                        count=count+1;
                    end
                % Inner Node
                elseif tempNodeType(j)==2
                    v1=zeros(1,3);
                    v2=zeros(1,3);
                    if j==1
                        v1=node0(C(N),:)-node0(C(1),:);
                        v2=node0(C(2),:)-node0(C(1),:);
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        crease1Num=SearchCreaseNum(oldCreaseConnect,C(N),C(1));
                        crease2Num=SearchCreaseNum(oldCreaseConnect,C(2),C(1));
                    elseif j==N
                        v1=node0(C(1),:)-node0(C(N),:);
                        v2=node0(C(N-1),:)-node0(C(N),:);
                        v1=v1/norm(v1);
                        v2=v2/norm(v2); 
                        crease1Num=SearchCreaseNum(oldCreaseConnect,C(1),C(N));
                        crease2Num=SearchCreaseNum(oldCreaseConnect,C(N-1),C(N));
                    else
                        v1=node0(C(j+1),:)-node0(C(j),:);
                        v2=node0(C(j-1),:)-node0(C(j),:);
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        crease1Num=SearchCreaseNum(oldCreaseConnect,C(j+1),C(j));
                        crease2Num=SearchCreaseNum(oldCreaseConnect,C(j-1),C(j));
                    end
                    if norm(cross(v1,v2))~=0
                        crease1W=creaseWidthMat(crease1Num);
                        crease2W=creaseWidthMat(crease2Num);
                        
                        v1=v1/norm(v1);
                        v2=v2/norm(v2);
                        
                        vector=crease1W/2*v2+crease2W/2*v1;
                        
                        newNode(count,:)=node0(C(j),:)+vector;
                        newNode2OldNode(count)=C(j);
                        count=count+1;
                    elseif norm(cross(v1,v2))==0
                        newNode(count,:)=node0(C(j),:)+vector;
                        newNode2OldNode(count)=C(j);
                        count=count+1;
                    end
                end
            end

            % Formulate the connection information vector 
            % Formulate the bar length vector
            PanelOrder=zeros(1,N);    
            for j=1:N
               PanelOrder(j)=count-N-1+j;
               newCrease2OldCrease(count-N-1+j)=tempCreaseNum(j);
               if j==1
                   barInfo(count-N,:)=[count-N,count-1,count-N+1,count-2,i];
                   barConnect(count-N,:)=[count-N,count-1];      
               elseif j==2
                   barInfo(count-N-1+j,:)=[count-N,count-N+1,count-1,count-N+2,i];
                   barConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
               elseif j==N
                   barInfo(count-N-1+j,:)=[count-N-2+j,count-N-1+j,count-N-3+j,count-N,i];  
                   barConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
               else
                   barInfo(count-N-1+j,:)=[count-N-2+j,count-N-1+j,count-N-3+j,count-N+j,i];
                   barConnect(count-N-1+j,:)=[count-N-2+j,count-N-1+j];
               end
            end

            % Assign Areas to bars of Panels (Peripheral)
            if N==2
            else        
                tempNode=zeros(1,3);
                for j=1:N
                    tempNode=tempNode+node0(C(j),:);
                end
                centerNodetemp=tempNode/N;  
            end  
            newPanel{i}=PanelOrder;    


        end
        barType=ones(count-1,1);
        type1BarNum=count-1;

        %% Generate the additional structure
        A=size(newCrease2OldCrease);
        NewCreaseNum=A(1);
        panelCount=newPanelNum+1;
        BarNum=count;

        for i=1:oldCreaseNum
            if oldCreaseType(i)==1
            else
                countTemp=0;
                edge1=0;
                edge2=0;
                for j=1:NewCreaseNum
                    if i==newCrease2OldCrease(j)
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

                Node11=newNode(barInfo(edge1,1),:);
                Node12=newNode(barInfo(edge1,2),:);
                Node21=newNode(barInfo(edge2,1),:);
                Node22=newNode(barInfo(edge2,2),:);

                Node11Ref1=newNode(barInfo(edge1,3),:)-newNode(barInfo(edge1,1),:);
                Node11Ref2=newNode(barInfo(edge1,2),:)-newNode(barInfo(edge1,1),:);

                Node21Ref1=newNode(barInfo(edge1,4),:)-newNode(barInfo(edge1,2),:);
                Node21Ref2=newNode(barInfo(edge1,1),:)-newNode(barInfo(edge1,2),:);

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

                % Generate Three new node and their location for the analysis
                if Flag2D3D==3
                    newNode(count,:)=Node11-RefVector1;
                    newNode2OldNode(count)=0;
                    count=count+1;
                    newNode(count,:)=Node12-RefVector2;
                    newNode2OldNode(count)=0;
                    count=count+1;
                    newNode(count,:)=0.5*(Node11+Node12)-0.5*(RefVector1+RefVector2);
                    newNode2OldNode(count)=0;
                    count=count+1;
                elseif Flag2D3D==2
                    newNode(count,:)=(Node11+Node21)/2;
                    newNode2OldNode(count)=0;
                    count=count+1;
                    newNode(count,:)=(Node12+Node22)/2;
                    newNode2OldNode(count)=0;
                    count=count+1;
                    newNode(count,:)=0.25*(Node12+Node22+Node11+Node21);
                    newNode2OldNode(count)=0;
                    count=count+1;
                end

                % Generate New Panels added to the structure
                newPanel{panelCount}=[count-1 count-3 NodeIndex11];
                newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                newPanel{panelCount}=[NodeIndex21 count-3 count-1];
                newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                newPanel{panelCount}=[NodeIndex21 count-1 NodeIndex22];
                newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                newPanel{panelCount}=[NodeIndex22 count-2 count-1];
                newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                newPanel{panelCount}=[NodeIndex12 count-1 count-2];
                newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;
                newPanel{panelCount}=[NodeIndex12 count-1 NodeIndex11];
                newPanel2OldPanel(panelCount)=0;
                panelCount=panelCount+1;

                barType(BarNum)=2;
                barConnect(BarNum,:)=[NodeIndex11,count-3];
                sprIJKL(BarNum,:)=[0,0,0,0]; 
                BarNum=BarNum+1;

                barType(BarNum)=3;
                barConnect(BarNum,:)=[NodeIndex11,count-1];
                sprIJKL(BarNum,:)=[NodeIndex12,count-1,NodeIndex11,count-3];
                newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                barType(BarNum)=3;
                barConnect(BarNum,:)=[NodeIndex12,count-1];
                sprIJKL(BarNum,:)=[count-2,count-1,NodeIndex12,NodeIndex11];
                newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                barType(BarNum)=2;
                barConnect(BarNum,:)=[NodeIndex12,count-2];
                sprIJKL(BarNum,:)=[0,0,0,0];            
                BarNum=BarNum+1;

                barType(BarNum)=4;
                barConnect(BarNum,:)=[count-3,count-1];
                sprIJKL(BarNum,:)=[NodeIndex11,count-1,count-3,NodeIndex21];  
                newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                barType(BarNum)=4;
                barConnect(BarNum,:)=[count-2,count-1];
                sprIJKL(BarNum,:)=[NodeIndex12,count-2,count-1,NodeIndex22]; 
                newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                barType(BarNum)=2;
                barConnect(BarNum,:)=[NodeIndex21,count-3];
                sprIJKL(BarNum,:)=[0 0 0 0];
                BarNum=BarNum+1;

                barType(BarNum)=3;
                barConnect(BarNum,:)=[NodeIndex21,count-1];
                sprIJKL(BarNum,:)=[count-3,count-1,NodeIndex21,NodeIndex22];
                newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                barType(BarNum)=3;
                barConnect(BarNum,:)=[NodeIndex22,count-1];
                sprIJKL(BarNum,:)=[count-2,NodeIndex22,count-1,NodeIndex21];
                newCrease2OldCrease(BarNum)=i;
                BarNum=BarNum+1;

                barType(BarNum)=2;
                barConnect(BarNum,:)=[NodeIndex22,count-2];
                sprIJKL(BarNum,:)=[0 0 0 0];
                BarNum=BarNum+1;  
            end 
        end

        %% Add bars for the bending of panels
        centerNodeStart=count-1;
        panelInnerBarStart=BarNum-1;
        for i=1:newPanelNum
            B=newPanel(i);
            C=cell2mat(B);
            D=size(C);
            N=D(2);

            if N==2
            else        
                tempNode=zeros(1,3);
                for j=1:N
                    tempNode=tempNode+newNode(C(j),:);
                end
                newNode(count,:)=tempNode/N;   
                newNode2OldNode(count)=-1;
                count=count+1;

                for j=1:N
                    barType(BarNum)=5;
                    barConnect(BarNum,:)=[C(j),count-1];
                    if j==1
                        sprIJKL(BarNum,:)=[C(N),C(1),count-1,C(2)];
                        newPanel{panelCount}=[count-1,C(1),C(2)];
                        newPanel2OldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    elseif j==N
                        sprIJKL(BarNum,:)=[C(N-1),C(N),count-1,C(1)];
                        newPanel{panelCount}=[count-1,C(1),C(N)];
                        newPanel2OldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    else
                        sprIJKL(BarNum,:)=[C(j-1),C(j),count-1,C(j+1)];
                        newPanel{panelCount}=[count-1,C(j),C(j+1)];
                        newPanel2OldPanel(panelCount)=i;
                        panelCount=panelCount+1;
                    end            
                    BarNum=BarNum+1;        
                end  
            end   
        end

        %% Calculate IJKL for panel bars
        count=1;
        for i=1:oldCreaseNum
            if oldCreaseType(i)==1
            else
                countTemp=0;
                edge1=0;
                edge2=0;
                for j=1:NewCreaseNum
                    if i==newCrease2OldCrease(j)
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

                if norm(newNode(Node11,:)-newNode(Node21,:))<norm(newNode(Node11,:)-newNode(Node22,:))
                else
                    tempNode=Node22;
                    Node22=Node21;
                    Node21=tempNode;
                end

                Panel1=barInfo(edge1,5);
                Panel2=barInfo(edge2,5);
                CenterNode=type1BarNum+3*(count-1)+3;

                % These are rotational springs at peripheral of panels, they use crease properties.
                sprIJKL(edge1,:)=[CenterNode,Node11,Node12,centerNodeStart+Panel1];
                sprIJKL(edge2,:)=[CenterNode,Node22,Node21,centerNodeStart+Panel2];
                count=count+1;
            end
        end
        
        
    %% If using the concentrated hinge model         
    else 

        A=size(node0);
        B=size(panel0);
        newPanelNum=B(2);
        nodeNum=A(1);
        centerNodeStart=nodeNum;
                
        newNode=zeros(newPanelNum+nodeNum,3);
        newNode(1:nodeNum,1:3)=node0(1:nodeNum,1:3);
        
        newNode2OldNode=zeros(newPanelNum+nodeNum,1);
        for i=1:nodeNum
            newNode2OldNode(i)=i;
        end
        
        panelCount=newPanelNum+1;
        newPanel=panel0;
        
        barCount=1;
        barConnect=zeros(5,2);
        barType=zeros(5,1);
        
        % This array stores how many time the bar appears, 
        % Value greater than 1 indicates The bar is shared by panels
        barApperenceNum=ones(5,1);
        
        sprIJKL=zeros(5,4);
        newCrease2OldCrease=zeros(5,1);  
        
        % Fromulation of Panel and Node
        for i=1:newPanelNum
            B=panel0(i);
            C=cell2mat(B);
            N=size(C,2);     
            centerNode=zeros(1,3);           
            
            for j=1:N
                % Form the coordinates of center node
                centerNode=centerNode+newNode(C(j),:);    
                % Form the panel information
                if j==1
                    newPanel{panelCount}=[nodeNum+i,C(1),C(N)];
                    newPanel2OldPanel(panelCount)=i;
                else
                    newPanel{panelCount}=[nodeNum+i,C(j),C(j-1)];
                    newPanel2OldPanel(panelCount)=i;
                end 
                panelCount=panelCount+1;
            end
            centerNode=centerNode/N;
            newNode(i+nodeNum,:)=centerNode;
        end
        
        % Calculate information for type 1 bars (at the edge)
        for i=1:newPanelNum
            B=panel0(i);
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
                       if (NodeMin==barConnect(k,1) && NodeMax==barConnect(k,2))
                           check=1;
                           barIndex=k;
                           barApperenceNum(barIndex)=barApperenceNum(barIndex)+1;
                       end
                    end
                    if check==0
                       barConnect(barCount,:)=[NodeMin,NodeMax];
                       barType(barCount)=1;
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
                       if (NodeMin==barConnect(k,1) && NodeMax==barConnect(k,2))
                           check=1;
                           barIndex=k;
                           barApperenceNum(barIndex)=barApperenceNum(barIndex)+1;
                       end
                    end
                    if check==0
                       barConnect(barCount,:)=[NodeMin,NodeMax];
                       barType(barCount)=1;
                       barCount=barCount+1;
                       barApperenceNum(barCount)=1;
                    else
                    end
                end
            end
        end
        
        panelInnerBarStart=barCount-1;        
        type1BarNum=barCount-1;   

        % Calculate information for type 5 bars (bars inside the panels)
        for i=1:newPanelNum
            B=panel0(i);
            C=cell2mat(B);
            N=size(C,2);             
            
            for j=1:N  
                barConnect(barCount,:)=[centerNodeStart+i,C(j)];               
                barType(barCount)=5;
                barApperenceNum(barCount)=1;       
                if j==1
                    sprIJKL(barCount,:)=[C(2),C(1),nodeNum+i,C(N)];
                elseif j==N
                    sprIJKL(barCount,:)=[C(1),C(N),nodeNum+i,C(N-1)];
                else
                    sprIJKL(barCount,:)=[C(j+1),C(j),nodeNum+i,C(j-1)];                    
                end
                barCount=barCount+1;
            end
        end     
        
        % Formulate IJKL matrix for bars connecting different panels
        for i=1:barCount-1
            if barApperenceNum(i)>1
                nodeMin=barConnect(i,1);
                nodeMax=barConnect(i,2);
                temp=1;
                CenterIndex=zeros(2,1);
                
                for j=1:newPanelNum
                     B=panel0(j);
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
                sprIJKL(i,:)=[nodeNum+CenterIndex(1),nodeMin,nodeMax,nodeNum+CenterIndex(2)];    
            end
        end
    end
    newPanel=newPanel(newPanelNum+1:end);
    newPanel2OldPanel=newPanel2OldPanel(newPanelNum+1:end);
    
    A=size(barConnect);
    barNum=A(1);
    newCrease2OldCrease(panelInnerBarStart:barNum)=zeros(size(panelInnerBarStart:barNum));
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