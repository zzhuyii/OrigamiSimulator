% This function solves the target rotaion angle of creases after heating
% using the resolved temperature profile of the origami

function targetRot=Thermal_UpdateRotZeroStrain(obj,T,thermal)

    targetRot=zeros(obj.oldCreaseNum,1);
    for i=1:obj.oldCreaseNum
        if obj.creaseRef(i,3)~=0
            
            bar1=obj.creaseRef(i,1);
            bar2=obj.creaseRef(i,2);        
            bar5=obj.creaseRef(i,5);
            bar6=obj.creaseRef(i,6);

            % obtain nodal temperature on the edge spring line
            node1=obj.barConnect(bar1,1);
            node2=obj.barConnect(bar1,2);        
            node3=obj.barConnect(bar2,1);
            node4=obj.barConnect(bar2,2);

            T1=T(node1)-obj.currentT(node1);   
            T2=T(node2)-obj.currentT(node2);
            T3=T(node3)-obj.currentT(node3);
            T4=T(node4)-obj.currentT(node4);

            % heating applied on the center spring line
            node5=obj.barConnect(bar5,1);
            node6=obj.barConnect(bar5,2);        
            node7=obj.barConnect(bar6,1);
            node8=obj.barConnect(bar6,2);
            
            T5=0;
            T6=0;
            T7=0;
            
            if node5==node6
                T5=T(node5)-obj.currentT(node5);
                T6=T(node7)-obj.currentT(node7);
                T7=T(node8)-obj.currentT(node8);
            elseif node5==node7
                T5=T(node5)-obj.currentT(node5);
                T6=T(node6)-obj.currentT(node6);
                T7=T(node8)-obj.currentT(node8);
            elseif node5==node8
                T5=T(node5)-obj.currentT(node5);
                T6=T(node7)-obj.currentT(node7);
                T7=T(node6)-obj.currentT(node6);
            elseif node6==node7
                T5=T(node6)-obj.currentT(node6);
                T6=T(node5)-obj.currentT(node5);
                T7=T(node8)-obj.currentT(node8);
            elseif node6==node8
                T5=T(node6)-obj.currentT(node6);
                T6=T(node5)-obj.currentT(node5);
                T7=T(node7)-obj.currentT(node7);
            elseif node7==node8
                T5=T(node7)-obj.currentT(node7);
                T6=T(node5)-obj.currentT(node5);
                T7=T(node6)-obj.currentT(node6);
            end
            
            dTave = 1/2*1/3*(T5+T6+T7) + 1/2*1/4*(T1+T2+T3+T4);
            deltaAlpha=thermal.deltaAlpha(i);
            rot = obj.Thermal_Timoshenko(dTave, deltaAlpha, thermal,...
                                         obj.creaseWidthVec(i));
            
            targetRot(i)=rot;
        end
    end  
end