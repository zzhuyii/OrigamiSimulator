%% Assemble stiffness of bars
% This function assembles the global stiffness matrix of the bars
%
% Input: 
%       U: current deformation field;
%       Sx: the second PK stress for bar elements;
%       C: tangent stiffness of bar elements;
%       barArea: area of bars;
%       barLength: length of bars;
%       barConnect: the two nodal number of bars;
%       newNode: nodal coordinates of meshed origami;
% Output: 
%       Kbar: global stiffness matrix from bars;
%

function [Kbar]=Bar_GlobalStiffAssemble(obj,U,Sx,C,barArea,barLength,barConnect,newNode)

    A=size(newNode);
    N_node=A(1); % Number of nodes
    
    A=size(C);
    N=A(1); % Number of bars
    
    Kbar=zeros(3*N_node,3*N_node);
%     Kbar_1=zeros(3*N_node,3*N_node);
%     Ktemp_his={};
%     B2_his={};
%     K_temp_part1_his={};
%     B1_B2_U_his={};

  
     %% This is the non-vectorized code
    
%     for i=1:N
%        NodeIndex1=barConnect(i,1);
%        NodeIndex2=barConnect(i,2);
%        node1=newNode(NodeIndex1,:);
%        node2=newNode(NodeIndex2,:);
% 
%        B1=1/(barLength(i)^2)*[-(node2-node1) (node2-node1)];
%        iden=eye(3);
%        B2=1/(barLength(i)^2)*[iden -iden; -iden iden];
% 
%        Utemp=[U(NodeIndex1,:)';U(NodeIndex2,:)'];   
%        Ktemp=C(i)*barArea(i)*barLength(i)*(B1'+B2*Utemp)*(B1'+B2*Utemp)'+ ...
%            Sx(i)*barArea(i)*barLength(i)*B2;
% 
%        index1=3*NodeIndex1-2;
%        index2=3*NodeIndex2-2;
% 
%        Kbar(index1:index1+2,index1:index1+2)=Kbar(index1:index1+2,index1:index1+2)+Ktemp(1:3,1:3);
%        Kbar(index2:index2+2,index2:index2+2)=Kbar(index2:index2+2,index2:index2+2)+Ktemp(4:6,4:6);
%        Kbar(index1:index1+2,index2:index2+2)=Kbar(index1:index1+2,index2:index2+2)+Ktemp(1:3,4:6);
%        Kbar(index2:index2+2,index1:index1+2)=Kbar(index2:index2+2,index1:index1+2)+Ktemp(4:6,1:3);
%     end
%     
    %% Vectorized code
    
    NodeIndex1_temp=barConnect(:,1);
    NodeIndex2_temp=barConnect(:,2);
    node1_temp=newNode(NodeIndex1_temp,:);
    node2_temp=newNode(NodeIndex2_temp,:);
    
    barLength_square=barLength.*barLength;
    
    B1_temp=(1./barLength_square).*([-(node2_temp-node1_temp), (node2_temp-node1_temp)]);
    iden=eye(3);
    B2_pattern=[iden -iden; -iden iden];
    
    B2_temp=[1./barLength_square.*B2_pattern(1,:),...
        1./barLength_square.*B2_pattern(2,:),...
        1./barLength_square.*B2_pattern(3,:),...
        1./barLength_square.*B2_pattern(4,:),...
        1./barLength_square.*B2_pattern(5,:),...
        1./barLength_square.*B2_pattern(6,:)];
        
    U_temp=[U(NodeIndex1_temp,:) U(NodeIndex2_temp,:)];
    
    B2_U=[dot(B2_temp(:,1:6),U_temp,2),...
          dot(B2_temp(:,7:12),U_temp,2),...
          dot(B2_temp(:,13:18),U_temp,2),...
          dot(B2_temp(:,19:24),U_temp,2),...
          dot(B2_temp(:,25:30),U_temp,2),...
          dot(B2_temp(:,31:36),U_temp,2)];
      
    B1_B2_U = B1_temp + B2_U;
    
    K_temp_part1=[B1_B2_U(:,1).*B1_B2_U,...
                  B1_B2_U(:,2).*B1_B2_U,...
                  B1_B2_U(:,3).*B1_B2_U,...
                  B1_B2_U(:,4).*B1_B2_U,...
                  B1_B2_U(:,5).*B1_B2_U,...
                  B1_B2_U(:,6).*B1_B2_U];
              
    K_temp=C.*barArea.*barLength.*K_temp_part1+Sx.*barArea.*barLength.*B2_temp;    
      
%     Ktemp=C(i)*barArea(i)*barLength(i)*(B1'+B2*Utemp)*(B1'+B2*Utemp)'+ ...
%            Sx(i)*barArea(i)*barLength(i)*B2;  

    index1=3*NodeIndex1_temp-2;
    index2=3*NodeIndex2_temp-2;

    index_dim=[index1, index1+1, index1+2,...
         index2, index2+1, index2+2];
    
    index=zeros(N,6,6,2);
    
    index(:,1,:,2)=index_dim;    
    index(:,2,:,2)=index_dim;
    index(:,3,:,2)=index_dim;
    index(:,4,:,2)=index_dim;
    index(:,5,:,2)=index_dim;
    index(:,6,:,2)=index_dim;
    
    index(:,:,1,1)=index_dim;    
    index(:,:,2,1)=index_dim;
    index(:,:,3,1)=index_dim;
    index(:,:,4,1)=index_dim;
    index(:,:,5,1)=index_dim;
    index(:,:,6,1)=index_dim;
    
    index=reshape(index, [N*36,2]);    
    K_temp=K_temp(:);
    for i=1:length(K_temp)
        Kbar(index(i,1),index(i,2))=Kbar(index(i,1),index(i,2))+K_temp(i);
    end
    
%     t=15;
%     reshape(K_temp(t,:),[6,6])
%     Ktemp_his{t}
%     
%     B1_B2_U(t,:)
%     B1_B2_U_his{t}
%     
%     reshape(K_temp_part1(t,:),[6,6])
%     K_temp_part1_his{t}
    
end