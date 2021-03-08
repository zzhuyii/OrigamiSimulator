%% Identify Creases
% This part of the code is used to indentify and number the creases based
% on the input data, so that the model is ready for meshing; 
%
% Input include:
%    {panel0}: stores the node number of each panel; 
%    [node0]: stores the coordinates of each nodes;
%
% Output include: 
%    [oldCreaseNum]: total number of creases;
%    [oldCreaseConnect]: two node numbers of the creases;
%    [oldCreaseType]: =1 for boundary crease; !=1 for middle crasese;
%
        

function Mesh_AnalyzeOriginalPattern(obj)

    A=size(obj.node0);
    B=size(obj.panel0);
    oldPanelNum=B(2);
    NodeNum=A(1);

    obj.oldCreaseNum=0;
    % Num of creases of the old pattern
    obj.oldCreaseConnect=zeros(5,2);
    % Stores the two nodes of a given crease
    obj.oldCreaseType=zeros(5,1);
    % Type of creases of the old pattern

    for i=1:oldPanelNum
        B=obj.panel0(i);
        C=cell2mat(B);
        D=size(C);
        N=D(2);

        % Identifying Creases of the original pattern
        % Start with identifying the first creases
        if obj.oldCreaseNum==0
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N));
            obj.oldCreaseConnect(1,:)=[minNode maxNode];
            obj.oldCreaseNum=obj.oldCreaseNum+1;
            obj.oldCreaseType(1)=1;        
        else
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N)); 
            A=[minNode maxNode];
            check=0;        
            for j=1:obj.oldCreaseNum
                if eq(A,obj.oldCreaseConnect(j,:))
                    check=check+1;
                    obj.oldCreaseType(j)=obj.oldCreaseType(j)+1;
                end
            end        
            if check ==0
               obj.oldCreaseNum=obj.oldCreaseNum+1;
               obj.oldCreaseConnect(obj.oldCreaseNum,:)=A;
               obj.oldCreaseType(obj.oldCreaseNum)=1;
            end        
        end
        % identifying remaining creases
        for j=1:N-1
            minNode=min(C(j),C(j+1));
            maxNode=max(C(j),C(j+1)); 
            A=[minNode maxNode];
            check=0;        
            for k=1:obj.oldCreaseNum
                if eq(A,obj.oldCreaseConnect(k,:))
                    check=check+1;
                    obj.oldCreaseType(k)=obj.oldCreaseType(k)+1;
                end
            end        
            if check ==0
               obj.oldCreaseNum=obj.oldCreaseNum+1;
               obj.oldCreaseConnect(obj.oldCreaseNum,:)=A;
               obj.oldCreaseType(obj.oldCreaseNum)=1;
            end 
        end
    end
    
    obj.creaseWidthVec=zeros(obj.oldCreaseNum,1);
    obj.creaseThickVec=zeros(obj.oldCreaseNum,1);
end



