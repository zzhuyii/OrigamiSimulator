%% Identify Creases
% This part of the code is used to indentify and number the creases based
% on the input data 
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
%

function[oldCreaseNum,oldCreaseConnect,oldCreaseType]=Mesh_IdentifyCrease(node0,panel0)

    A=size(node0);
    B=size(panel0);
    oldPanelNum=B(2);
    NodeNum=A(1);

    oldCreaseNum=0;
    % Num of creases of the old pattern
    oldCreaseConnect=zeros(5,2);
    % Stores the two nodes of a given crease
    oldCreaseType=zeros(5,1);
    % Type of creases of the old pattern

    for i=1:oldPanelNum
        B=panel0(i);
        C=cell2mat(B);
        D=size(C);
        N=D(2);

        % Identifying Creases of the original pattern
        % Start with identifying the first creases
        if oldCreaseNum==0
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N));
            oldCreaseConnect(1,:)=[minNode maxNode];
            oldCreaseNum=oldCreaseNum+1;
            oldCreaseType(1)=1;        
        else
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N)); 
            A=[minNode maxNode];
            check=0;        
            for j=1:oldCreaseNum
                if eq(A,oldCreaseConnect(j,:))
                    check=check+1;
                    oldCreaseType(j)=oldCreaseType(j)+1;
                end
            end        
            if check ==0
               oldCreaseNum=oldCreaseNum+1;
               oldCreaseConnect(oldCreaseNum,:)=A;
               oldCreaseType(oldCreaseNum)=1;
            end        
        end
        % identifying remaining creases
        for j=1:N-1
            minNode=min(C(j),C(j+1));
            maxNode=max(C(j),C(j+1)); 
            A=[minNode maxNode];
            check=0;        
            for k=1:oldCreaseNum
                if eq(A,oldCreaseConnect(k,:))
                    check=check+1;
                    oldCreaseType(k)=oldCreaseType(k)+1;
                end
            end        
            if check ==0
               oldCreaseNum=oldCreaseNum+1;
               oldCreaseConnect(oldCreaseNum,:)=A;
               oldCreaseType(oldCreaseNum)=1;
            end 
        end
    end
end



