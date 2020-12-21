% this function converts the crease heating to nodal heating

function Qload=Thermal_ConvertCreaseHeat2NodeHeat(obj,creaseHeating)
    nodeNum=size(obj.newNode);
    nodeNum=nodeNum(1);
    
    creaseNum=size(creaseHeating);
    creaseNum=creaseNum(1);
    
    Qload=zeros(obj.envLayer*nodeNum,1);
    
    for i=1:creaseNum
        bar1=obj.creaseRef(creaseHeating(i,1),1);
        bar2=obj.creaseRef(creaseHeating(i,1),2);        
        bar5=obj.creaseRef(creaseHeating(i,1),5);
        bar6=obj.creaseRef(creaseHeating(i,1),6);
        
        % heating applied on the side spring line
        node1=obj.barConnect(bar1,1);
        node2=obj.barConnect(bar1,2);        
        node3=obj.barConnect(bar2,1);
        node4=obj.barConnect(bar2,2);
                
        Qload(node1)=1/4*1/2*creaseHeating(i,2);
        Qload(node2)=1/4*1/2*creaseHeating(i,2);
        Qload(node3)=1/4*1/2*creaseHeating(i,2);
        Qload(node4)=1/4*1/2*creaseHeating(i,2);        
        
        % heating applied on the center spring line
        node5=obj.barConnect(bar5,1);
        node6=obj.barConnect(bar5,2);        
        node7=obj.barConnect(bar6,1);
        node8=obj.barConnect(bar6,2);
        
        if node5==node6
            Qload(node5)=1/2*2/3*creaseHeating(i,2);
            Qload(node7)=1/2*1/6*creaseHeating(i,2);
            Qload(node8)=1/2*1/6*creaseHeating(i,2);
        elseif node5==node7
            Qload(node5)=1/2*2/3*creaseHeating(i,2);
            Qload(node6)=1/2*1/6*creaseHeating(i,2);
            Qload(node8)=1/2*1/6*creaseHeating(i,2);
        elseif node5==node8
            Qload(node5)=1/2*2/3*creaseHeating(i,2);
            Qload(node7)=1/2*1/6*creaseHeating(i,2);
            Qload(node6)=1/2*1/6*creaseHeating(i,2);
        elseif node6==node7
            Qload(node6)=1/2*2/3*creaseHeating(i,2);
            Qload(node5)=1/2*1/6*creaseHeating(i,2);
            Qload(node8)=1/2*1/6*creaseHeating(i,2);
        elseif node6==node8
            Qload(node6)=1/2*2/3*creaseHeating(i,2);
            Qload(node5)=1/2*1/6*creaseHeating(i,2);
            Qload(node7)=1/2*1/6*creaseHeating(i,2);
        elseif node7==node8
            Qload(node7)=1/2*2/3*creaseHeating(i,2);
            Qload(node5)=1/2*1/6*creaseHeating(i,2);
            Qload(node6)=1/2*1/6*creaseHeating(i,2);
        end
    end    
end


