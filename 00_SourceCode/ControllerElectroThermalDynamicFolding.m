% This class stores the properties used to control the thermally/ electro
% thermally acivated loading of an origami

classdef ControllerElectroThermalDynamicFolding < handle
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
        
        % time increment of each step
        dt

        % total step
        step

        % matrix to store the target heating history;
        targetCreaseHeatingHis=[]
        
        % show the figure of a final deformed shape
        plotOpen=1

        % crop frames for video
        videoCropRate=100
        
        % plot animation
        videoOpen=1
        
        % plot details
        detailFigOpen=0    

        % Rayleigh damping
        alpha=0.0001
        beta=0.0001

        
        %% Thmeral boundary
        
        % Define the BC for thermal conduction, to confine the maximum air
        % thickness. The vector stores the number of panels that serves as the BC
        % for heat transfer (ie. those that are Si wafers).
        thermalBoundaryPanelVec=[]

        % This vector stores extra nodes that should be at RT to adjust the BC.
        roomTempNode=[]
        
        
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
        
        % The history of nodal forces
        FnodalHis
        
        % The history of bar strain
        barExHis
        
        % The history of bar stress
        barSxHis
        
        % The history of the rotational spring moment
        sprMHis
        
        % The history of the rotational spring rotation
        sprRotHis
    
    end    
end