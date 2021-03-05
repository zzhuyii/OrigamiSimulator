%% formulate the mass matrix

function NodalMass=Dynamic_MassVector(obj)

    NumNode=size(obj.newNode,1);
    NodalMass=zeros(NumNode,1);
    NumBar=size(obj.barConnect);
    
    gammaCrease=obj.densityCrease;
    gammaPanel=obj.densityPanel;
    
    NumPanelEnd=size(obj.newPanel,2);
    NumPanelStart=size(obj.panel0,2);
   
    for i=NumPanelStart+1:NumPanelEnd
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
       
       % Distribute the mass 
       if barType1==5 || barType2==5 || barType3==5 % panel triangle
           PanelThick=obj.panelThickMat(obj.newPanel2OldPanel(i));
           triangleMass=area*gammaPanel*PanelThick;
           if barType1==5 && barType2==5
               NodalMass(nodeVec(2))=...
                   1/3*triangleMass+NodalMass(nodeVec(2));
           else
               NodalMass(nodeVec(2))=...
                   1/3*triangleMass+NodalMass(nodeVec(2));
           end
           if barType2==5 && barType3==5
               NodalMass(nodeVec(3))=...
                   1/3*triangleMass+NodalMass(nodeVec(3));
           else
               NodalMass(nodeVec(3))=...
                   1/3*triangleMass+NodalMass(nodeVec(3));
           end 
           if barType1==5 && barType3==5
               NodalMass(nodeVec(1))=...
                   1/3*triangleMass+NodalMass(nodeVec(1));
           else
               NodalMass(nodeVec(1))=...
                   1/3*triangleMass+NodalMass(nodeVec(1));
           end 
       else % crease triangle
           oldCreaseNum=max([obj.newCrease2OldCrease(NumBar1),...
               obj.newCrease2OldCrease(NumBar2),...
               obj.newCrease2OldCrease(NumBar3)]);
           
           CreaseThick=obj.creaseThickMat(oldCreaseNum);
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