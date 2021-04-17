% This class stores the properties used to control the Newton-Raphson
% loading of an origami

classdef ControllerNRLoading  < handle
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
        
        % the applied load
        load
        % a sampel code is the following
        % loadForce=3;
        % load=[29,0,0,loadForce;
        %       30,0,0,loadForce;
        %       31,0,0,loadForce;
        %       32,0,0,loadForce;];
        
        % the total number of incremental steps
        increStep=50
        
        % the tolerance for each iteration
        tol=1*10^-5
        
        % the maximum allowed iteration number
        iterMax=30        
        
        % show the figure of a final deformed shape
        plotOpen=1
        
        % plot animation
        videoOpen=1
        
        % plot details
        detailFigOpen=0
        
        %% Output properties
        % the history of displacement field
        Uhis
        
        % the history of strain energy
        strainEnergyHis
        
        % the history of loading
        loadHis
        
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