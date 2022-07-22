%% formulate the mass matrix

function NodalMass=Dynamic_MassVector(obj)

    NumNode=size(obj.newNode,1);
    NodalMass=zeros(NumNode,1);
    NumBar=size(obj.barConnect);
    
    gammaCrease=obj.densityCrease;
    gammaPanel=obj.densityPanel;
    
    NumPanelEnd=size(obj.newPanel,2);
   
    for i=1:NumPanelEnd
       nodeVec=obj.newPanel{i};
       
       vec1=obj.newNode(nodeVec(2),:)-obj.newNode(nodeVec(1),:);
       vec2=obj.newNode(nodeVec(3),:)-obj.newNode(nodeVec(1),:);       
       area=0.5*norm(cross(vec1,vec2));
       
       node1bar1=min(nodeVec(1),nodeVec(2));
       node2bar1=max(nodeVec(1),nodeVec(2));
       NumBar1=1;
       for j=1:NumBar
           if obj.barConnect(j,1)==node1bar1 && obj.barConnect(j,2)==node2bar1
               NumBar1=j;
           end
       end
       barType1=obj.barType(NumBar1);
       
       node1bar2=min(nodeVec(2),nodeVec(3));
       node2bar2=max(nodeVec(2),nodeVec(3));
       NumBar2=1;
       for j=1:NumBar
           if obj.barConnect(j,1)==node1bar2 && obj.barConnect(j,2)==node2bar2
               NumBar2=j;
           end
       end
       barType2=obj.barType(NumBar2);
       
       node1bar3=min(nodeVec(1),nodeVec(3));
       node2bar3=max(nodeVec(1),nodeVec(3));
       NumBar3=1;
       for j=1:NumBar
           if obj.barConnect(j,1)==node1bar3 && obj.barConnect(j,2)==node2bar3
               NumBar3=j;
           end
       end
       barType3=obj.barType(NumBar3);
       
       % set up mass for non-compliant crease
       if obj.compliantCreaseOpen==0
           barType1=5;
           barType2=5;
           barType3=5;
       end
       
       % Distribute the mass 
       if barType1==5 || barType2==5 || barType3==5 % panel triangle
           PanelThick=obj.panelThickVec(obj.newPanel2OldPanel(i));
           triangleMass=area*gammaPanel*PanelThick;
           
           NodalMass(nodeVec(1))=...
               1/3*triangleMass+NodalMass(nodeVec(1));
           NodalMass(nodeVec(2))=...
               1/3*triangleMass+NodalMass(nodeVec(2));
           NodalMass(nodeVec(3))=...
               1/3*triangleMass+NodalMass(nodeVec(3));
           
       else % crease triangle
           if obj.compliantCreaseOpen==1
               oldCreaseNum=max([obj.newCrease2OldCrease(NumBar1),...
                   obj.newCrease2OldCrease(NumBar2),...
                   obj.newCrease2OldCrease(NumBar3)]);

               CreaseThick=obj.creaseThickVec(oldCreaseNum);
              
               triangleMass=area*gammaCrease*CreaseThick;

               NodalMass(nodeVec(1))=...
                   1/3*triangleMass+NodalMass(nodeVec(1));
               NodalMass(nodeVec(2))=...
                   1/3*triangleMass+NodalMass(nodeVec(2));
               NodalMass(nodeVec(3))=...
                   1/3*triangleMass+NodalMass(nodeVec(3));
           end
       end
    end
end