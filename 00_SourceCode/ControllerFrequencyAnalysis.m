% This class stores the properties used to control the Newton-Raphson
% loading of an origami

classdef ControllerFrequencyAnalysis  < handle
    properties
        % storing the support information
        supp
        
        % NonRigidSupport
        % 0 means the non rigid support is not activated.
        % 1 means the non rigid support is activated.
        nonRigidSupport=0        
        
        % first column stores node number
        % second column stores direction
        % third column stores stiffness
        suppElastic
        % a sample input here
        % suppElastic=[1,3,10000;4,3,10000]; 
        
        % show the figure of a final deformed shape
        plotOpen=1
        

        
    end
end