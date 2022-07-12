%% Assemble inner force of bars
% This function assembles the global inner force vector of bars 
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
%       Tbar: global inner force vector from bars;
%

function [Tbar]=Bar_GlobalForce(obj,U,Sx,C,barArea,barLength,barConnect,newNode)
    A=size(newNode);
    N=A(1);
    Tbar=zeros(3*N,1);    

    
%% Non-vectorized code 
%     A=size(C);
%     N=A(1);
%  
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
%        
%        Ttemp=Sx(i)*barArea(i)*barLength(i)*(B1'+B2*Utemp);
%        Ttemprecord{i}=Ttemp;
%        index1=3*NodeIndex1-2;
%        index2=3*NodeIndex2-2;
% 
%        Tbar(index1:index1+2)=Tbar(index1:index1+2)+Ttemp(1:3);
%        Tbar(index2:index2+2)=Tbar(index2:index2+2)+Ttemp(4:6);   
%     end
% 

%% Vectorized code

    NodeIndex1=barConnect(:,1);
    NodeIndex2=barConnect(:,2);
    node1=newNode(NodeIndex1,:);
    node2=newNode(NodeIndex2,:);
    
    B1n=(1./(barLength.*barLength)).*[-(node2-node1) (node2-node1)];  
    
    iden=eye(3);
    idenMat=[iden -iden; -iden iden];
    
    Utemp1=[U(NodeIndex1,:)';U(NodeIndex2,:)'];     
    B2Utemp=(1./(barLength.*barLength)).*(idenMat*Utemp1)';
    % Vectroized operation
    
    index1=3*NodeIndex1-2;
    index2=3*NodeIndex2-2;      

    Ttemp=Sx.*barArea.*barLength;
    Ttemp=Ttemp.*(B1n+B2Utemp);
    
    index=[index1,index1+1,index1+2,index2,index2+1,index2+2];
    index=index(:);
    Ttemp=Ttemp(:);    
    
    for i=1:length(index)
        Tbar(index(i))=Tbar(index(i))+Ttemp(i);
    end
    
end