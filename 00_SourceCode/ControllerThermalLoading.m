% This class stores the properties used to control the thermally/ electro
% thermally acivated loading of an origami

classdef ControllerThermalLoading < handle
    properties
        %% Input Proerperties
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
        
        % thermal loading step
        thermalStep
        
        % matrix to store the target heating. This is the incremental term 
        % of this step. First index stores the crease number the second
        % index stores the applied heating power;
        targetCreaseHeating      
        
        % the tolerance for each iteration
        tol
        
        % the maximum allowed iteration number
        iterMax=30
        
        % the total number of incremental steps for mechanically update the
        % equilibrium within each thermal increment
        increStep=1
        
        % show the figure of a final deformed shape
        plotOpen=1
        
        % plot animation
        videoOpen=1
        
        % plot details
        detailFigOpen=0       

        
        %% Thmeral boundary
        
        % Define the BC for thermal conduction, to confine the maximum air
        % thickness. The vector stores the number of panels that serves as the BC
        % for heat transfer (ie. those that are Si wafers).
        thermalBoundaryPanelMat

        % This vector stores extra nodes that should be at RT to adjust the BC.
        roomTempNode
        
        
        %% Timoshenko bimoph model constants
        % differential thermal coefficient
        deltaAlpha
        
        % E of mateiral 1
        Emat1
        
        % E of material 2
        Emat2
        
        % thickness of material 1
        tmat1
        
        % thickness of material 2
        tmat2
        
        
        %% Output properties
        % the history of displacement field
        Uhis
        
        % the history of strain energy
        strainEnergyHis
       
        % the history of temperature
        temperatureHis
    
    end    
end