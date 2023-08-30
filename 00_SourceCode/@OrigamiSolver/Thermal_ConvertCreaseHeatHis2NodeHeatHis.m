% this function converts the crease heating history to nodal heating
% history

function QloadHis=Thermal_ConvertCreaseHeatHis2NodeHeatHis(obj,creaseHeatingHis)
    nodeNum=size(obj.newNode);
    nodeNum=nodeNum(1);
    
    creaseNum=size(creaseHeatingHis);
    stepSize=creaseNum(2)-1;
    creaseNum=creaseNum(1);
    
    QloadHis=zeros((obj.envLayer+1)*nodeNum,stepSize);
    
    for i=1:creaseNum
        bar1=obj.creaseRef(creaseHeatingHis(i,1),1);
        bar2=obj.creaseRef(creaseHeatingHis(i,1),2);        
        bar5=obj.creaseRef(creaseHeatingHis(i,1),5);
        bar6=obj.creaseRef(creaseHeatingHis(i,1),6);
        
        % heating applied on the side spring line
        node1=obj.barConnect(bar1,1);
        node2=obj.barConnect(bar1,2);        
        node3=obj.barConnect(bar2,1);
        node4=obj.barConnect(bar2,2);
                
        QloadHis(node1,:)=QloadHis(node1)+1/4*1/2*creaseHeatingHis(i,2:stepSize+1);
        QloadHis(node2,:)=QloadHis(node2)+1/4*1/2*creaseHeatingHis(i,2:stepSize+1);
        QloadHis(node3,:)=QloadHis(node3)+1/4*1/2*creaseHeatingHis(i,2:stepSize+1);
        QloadHis(node4,:)=QloadHis(node4)+1/4*1/2*creaseHeatingHis(i,2:stepSize+1);        
        
        % heating applied on the center spring line
        node5=obj.barConnect(bar5,1);
        node6=obj.barConnect(bar5,2);        
        node7=obj.barConnect(bar6,1);
        node8=obj.barConnect(bar6,2);
        
        if node5==node6
            QloadHis(node5,:)=1/2*2/3*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node7,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node8,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
        elseif node5==node7
            QloadHis(node5,:)=1/2*2/3*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node6,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node8,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
        elseif node5==node8
            QloadHis(node5,:)=1/2*2/3*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node7,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node6,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
        elseif node6==node7
            QloadHis(node6,:)=1/2*2/3*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node5,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node8,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
        elseif node6==node8
            QloadHis(node6,:)=1/2*2/3*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node5,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node7,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
        elseif node7==node8
            QloadHis(node7,:)=1/2*2/3*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node5,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
            QloadHis(node6,:)=1/2*1/6*creaseHeatingHis(i,2:stepSize+1);
        end
    end    
end


