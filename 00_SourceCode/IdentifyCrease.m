%% Identify Creases
% This part of the code is used to indentify and number the creases based
% on the input data file [Panel] and [Node]. Output file include CreaseNum,
% which indicate the total number of creases; Crease[:,2] Arrary, which
% give the two node numbers of the two creases; and CreaseType[:] which
% keep track of how many times a crease appear. 

% If CreaseType(i)==1, it is a boundary crease. If CreaseType(i)!=1, it is
% a crease in the middel of tessellation.

function[OldCreaseNum,OldCreaseConnect,OldCreaseType]=IdentifyCrease(Node0,Panel0)

    A=size(Node0);
    B=size(Panel0);
    PanelNum=B(2);
    NodeNum=A(1);

    OldCreaseNum=0;
    % Num of creases of the old pattern
    OldCreaseConnect=zeros(10,2);
    % Stores the two nodes of a given crease
    OldCreaseType=zeros(10,1);
    % Type of creases of the old pattern

    for i=1:PanelNum
        B=Panel0(i);
        C=cell2mat(B);
        D=size(C);
        N=D(2);

        % Identifying Creases of the original pattern
        % Start with identifying the first creases
        if OldCreaseNum==0
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N));
            OldCreaseConnect(1,:)=[minNode maxNode];
            OldCreaseNum=OldCreaseNum+1;
            OldCreaseType(1)=1;        
        else
            minNode=min(C(1),C(N));
            maxNode=max(C(1),C(N)); 
            A=[minNode maxNode];
            check=0;        
            for j=1:OldCreaseNum
                if eq(A,OldCreaseConnect(j,:))
                    check=check+1;
                    OldCreaseType(j)=OldCreaseType(j)+1;
                end
            end        
            if check ==0
               OldCreaseNum=OldCreaseNum+1;
               OldCreaseConnect(OldCreaseNum,:)=A;
               OldCreaseType(OldCreaseNum)=1;
            end        
        end
        % identifying remaining creases
        for j=1:N-1
            minNode=min(C(j),C(j+1));
            maxNode=max(C(j),C(j+1)); 
            A=[minNode maxNode];
            check=0;        
            for k=1:OldCreaseNum
                if eq(A,OldCreaseConnect(k,:))
                    check=check+1;
                    OldCreaseType(k)=OldCreaseType(k)+1;
                end
            end        
            if check ==0
               OldCreaseNum=OldCreaseNum+1;
               OldCreaseConnect(OldCreaseNum,:)=A;
               OldCreaseType(OldCreaseNum)=1;
            end 
        end
    end
end



