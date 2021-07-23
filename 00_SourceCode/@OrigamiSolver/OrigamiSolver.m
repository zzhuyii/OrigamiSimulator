%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Yi Zhu, and Evgueni T. Filipov
%
% Discription: This code package implement a bar and hinge model based 
% simulator for active origami structures. The code package can capture 
% the following behaviors of active origami systems;
%
% (1) Simulate panel contact in origami;
% (2) Simulate electro-thermal actuation for folding origami;
% (3) Provide 3 different solver for large deformation loading;
% (4) Provide elastic support;
% (5) Provide compliant crease bar and hinge model for origami;
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% origami simulator, the origami contact, compliant crease, electro-thermal
% model presented in this package. 
%
% Reference:
% [1] Y. Zhu, E. T. Filipov (2020). 'Rapid Multi-Physic Simulation for 
%     Electro-Thermal Origami Robotic Systems' (submitted)
% [2] Y. Zhu, E. T. Filipov (2020). 'A Bar and Hinge Model for Simulating 
%     Bistability in Origami Structures with Compliant Creases' Journal of 
%     Mechanisms and Robotics, 021110-1. 
% [3] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' Proc. R. Soc. A, 475: 
%     20190366.       
% [4] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [5] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [6] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' Proc. R. Soc. A 473: 
%     20170348. 
% [7] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%

% This gives the origami solver class, the main class of the package for
% solving origami folding under different loading situations. Different
% properties and methods are introduced briefly in this file. Both the
% properties and methods are divided into "user use" and "internal", where
% the first categories are those used by common users, while the internal
% parts are for developers. 

classdef OrigamiSolver < handle
    properties

%%
%%%%%%%%%%%%%%%%%%%%%%%  User Defined Variables  %%%%%%%%%%%%%%%%%%%%%%%%

        %% Origami Geometrical Properties        
        % nodal coordinates of the original origami
        node0; 
        
        % panel connectivity of the orignal origami
        panel0;
        
        % flag2D3D is used to determine how the additional crease structure
        % is generated,2D means center bars are genertaed through using an 
        % averaged vector, 3D means the center bars are located at the 
        % original positoin.
        % 3 3D, 2 2D
        flag2D3D=2;
        
        % compliantCreaseOpen is used to determine if compliant crease
        % model should be used or not
        % 1 means include compliant crease model
        % 0 means using concentrated hinge model        
        compliantCreaseOpen=1;
        
        % crease width to generate the geometry of the compliant 
        % crease origami;
        creaseWidthVec;
        
        % nodal coordinates of the meshed origami
        % sometime user is expected to directly edit it to correct the
        % minor imperfection from the automated meshing
        newNode
        
        
        %% Origami Mechanical Properties
        % Young's modulus of panel
        panelE
        
        % Young's modulus of creases
        creaseE
        
        % Poisson ratio of panel
        panelPoisson
        
        % Poisson ratio of crease
        creasePoisson
        
        % stores the thickness of panel
        panelThickVec
        
        % stroes the thickness of creases
        creaseThickVec
        
        % Diagonal Rate is the factor that determine how much diagonal 
        % springs are stiffer than horizontal ones.
        diagonalRate=4
        
        % assumed width for panel to calculate panel bending stiffness
        panelW
        
        % distributing target zero-stress folding between three 
        % rotational srping lines
        zeroStrianDistributionFactor=0.5       

        
        
        %% Origami Thermal Properties       
        % Thermal conductivity vector for panel
        panelThermalConductVec
        
        % Thermal conductivity of crease
        creaseThermalConduct
        
        % Thermal conductivity ofsubmerged environment
        envThermalConduct
        
        % thickness of the submerged environment at RT
        t2RT
        
        % number of evironment layer
        envLayer=10;
        
        % Thermal dissipation angle
        thermalDissipation=20/180*3.14;
        
        % temperature of the submerged environment
        RT=21;
                
        
        %% Panel Contact Properties        
        % contact open
        % 1: means consider panel contact;
        % 0: means ignore panel contact
        contactOpen=0
        
        % ke: used to scale the magnitude of potentil
        ke

        % d0edge: d0 for points at the edge
        d0edge

        % d0center: d0 for points at the center
        d0center
        
        
        %% Dynamic Properties
        % density of crease
        densityCrease
        
        % density of panel
        densityPanel               
        
        
        %% Ploting properties        
        % the viewing angle for plotting
        viewAngle1=15;
        viewAngle2=30;
        
        % the range of the axis
        displayRange=1;
        
        % the ratio of the displayed range in the negative region
        displayRangeRatio=0.5;
        
        % whether to show the numbers
        showNumber = 1;
        
        % the color of faces for showing the numbering
        faceColorNumbering = [0.5,0.1,0.1];
        
        % the color of faces when ploting animation and deformed
        % configuration
        faceColorAnimation = 'yellow';
        
        % transparency of faces for showing the numbering
        faceAlphaNumbering = 0;
        
        % transparency of faces when ploting animation and deformed
        % configuration
        faceAlphaAnimation = 1;
        
        %% Origami Loading Controller
        loadingController={}
        
        % This controls if the solver will continue the solving to the
        % next sequence without initialization
        continuingLoading=0;
        
%%
%%%%%%%%%%%%%%%%%%%%%%%% Internal Varialbles %%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %% Origami Geometrical Properties     
        % the total number of creases in the original origami
        oldCreaseNum;
        
        % the two node of each crease in the original origami
        oldCreaseConnect;
        
        % the type of each crease in the original origami
        oldCreaseType;
        
        % the panel connectivity of the meahsed origami
        newPanel
        
        % [barType]: the type of the new bars;
        %       [1] Panel Bars
        %       [2] Vertical Crease Bars that connect different panels
        %       [3] Diagonal Crease Bars
        %       [4] Horizontal Crease Bars
        %       [5] Panel Diagonal Bars
        barType
        
        % [barConnect]: the two nodes that connects a bar;
        barConnect
        
        % conentivity information of the rotaional springs;
        sprIJKL        
        
        % type1BarNum: total number of panel bars;
        type1BarNum
        
        % panelInerBarStart: Number of the first type [5] bars;
        panelInnerBarStart
        
        % centerNodeStart: Number of the first node at the center of paenl;
        centerNodeStart
        
        % [newNode2OldNode]: mapping between the new node and the old node;
        newNode2OldNode
        
        % [newCrease2OldCrease]: mapping between the new crease and the 
        % old crease;
        newCrease2OldCrease
   
        % [newPanel2OldPanel]: mapping between the new panel and the 
        % old panel;
        newPanel2OldPanel
        
        % panelNum: total number of new panels;
        newPanelNum
        
        % how ith old crease are further divided into different new 
        % rotaional springs;
        creaseRef
        
        % the length of each bar
        barLength
        
        
        %% Internal Mechanical Properties
        % the area of each bar
        barArea
        
        % the rotational stiffness of each spring
        sprK        
            
        % stress free angle of creases
        currentRotZeroStrain
        
        % the zero strain folding angle of a crease
        currentSprZeroStrain  
        
        % current displacement filed
        currentU        
       
        % current external loading applied on structre
        currentAppliedForce
        
        
        %% Internal thermal properties
        % current nodal temperature
        currentT
        
        % current heating energy input
        currentQ
        
    end
    
    methods
        
%%
%%%%%%%%%%%%%%%%%%%%%% Methods for User to call  %%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Meshing related fucntions
        
        % analyze the input node0 and panel0 for meshing
        Mesh_AnalyzeOriginalPattern(obj)    
        
        % Generate the meshed geometry of the system
        Mesh_Mesh(obj)
        
        
        %% Ploting related methods
        
        % Plot the un-meshed origami pattern
        Plot_UnmeshedOrigami(obj)
        
        % Plot the meshed origami pattern
        Plot_MeshedOrigami(obj)
        
        % Plot the deformed shape and deformation history of origami
        Plot_DeformedShape(obj,undeformedNode,deformNode)
        Plot_DeformedShapeTemp(obj,thermal,newNode,deformNode,T)
        Plot_DeformedHis(obj,undeformedNode,UhisLoading)
        Plot_DeformedHisTemp(obj,beforeLoadingNode,UhisThermal,temperatureHistory)
        
        % Plot the detail loading information
        Plot_LoadHis(obj,loadHis,UhisLoading)
        Plot_Energy(obj,UhisLoading,strainEnergyLoading)
        Plot_ContactForce(obj,undeformedNode,deformNode,contactForce)
        Plot_LoadAndReaction(obj,undeformedNode,deformNode,...
            load,supp,nodeForce,loadForce)
        
        
        %% Solve the loading process
        Solver_Solve(obj)
        
        
        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Internal Methods  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Meshing methods
        
        % generate creaseRef matrix used for calculating the contact
        Mesh_NumberingForContact(obj);

        % calcualte the barLength 
        Mesh_BarLength(obj);
        
        % Calculate the mechanial properties
        Mesh_MechanicalProperties(obj);
        
        % calculate the zero strain folding of rotational springs
        sprZeroStrain=Mesh_CalculateZeroStrainFolding(obj,rotZeroStrain);
        
        
        %% Bar elements
        [Tbar]=Bar_GlobalForce(obj,U,Sx,C,barArea,barLength,barConnect,newNode);
        [Kbar]=Bar_GlobalStiffAssemble(obj,U,Sx,C,barArea,barLength,barConnect,newNode);
        [Ex]=Bar_Strain (obj,U,newNode,barArea,barConnect,barLength);
        [Sx,C]= Bar_Cons(obj,barType,Ex,panelE,creaseE);
        
        
        %% Spr elements
        [M,sprKadj]= Spr_Cons(obj,sprTargeZeroStrain,theta,...
            sprK,creaseRef,oldCreaseNum,panelInnerBarStart,sprIJKL,...
            newNode,U,compliantCreaseOpen);
        [Tspr]=Spr_GlobalForce(obj,U,M,sprIJKL,sprKadj,newNode);
        [Kspr]=Spr_GlobalStiffAssemble(obj,U,M,sprIJKL,sprKadj,newNode);
        [theta]=Spr_Theta(obj,U,sprIJKL,newNode);
        
        
        %% Contact potential
        [Tcontact,Kcontact]=Contact_AssembleForceStiffness(obj,...
            panel0,newNode2OldNode,newNode,U,ke,d0edge,d0center,...
            centerNodeStart,compliantCreaseOpen);
        [Dd2x2,Ddx]=Contact_DerivativeZone0(obj,Point,T1,T2,T3);
        [Dd2x2,Ddx]=Contact_DerivativeZone1(obj,Point,T1,d);
        [Dd2x2,Ddx]=Contact_DerivativeZone2(obj,Point,T1,T2);
        [d,zone,N1,N2]=Contact_P2TDistance(obj,Point,T1,T2,T3);
        
        
        %% Thermal models        
        % this function converts the crease heating to nodal heating
        Qload=Thermal_ConvertCreaseHeat2NodeHeat(obj,creaseHeating);
        
        % method to update rotation zero strain for thermal loading
        targetRot=Thermal_UpdateRotZeroStrain(obj,T,thermal);
        
        % to calculate the folding angle using Timoshenko Model
        rot=Thermal_Timoshenko(obj,Tcrease,deltaAlpha,thermal,creaseWidth);
        
        % assemble the heat conductivity matrix
        [thermalMat]=Thermal_AssembleConductMat(obj,thermal,U);
        
        % solve the temperature profile with heating energy
        [T,indexArray]=Thermal_SolveTemperature(obj,qin,thermalMat,thermalNodeNum,thermal);
        
        
        %% solver for loding               
        % solver for thermal loading
        [U,UhisThermal,energyHisThermal,temperatureHistory,...
            rotTargetZeroStrain,sprTargetZeroStrain]...
            =Solver_LoadingThermal(obj,thermal);     
        
        % Solver for NR loading
        [U,UhisLoading,loadHis,strainEnergyLoading,...
            nodeForce,loadForce,contactForce]=Solver_LoadingNR(obj,nr);
        
        % Solver for MGDCM loading
        [U,UhisLoading,loadHis,strainEnergyLoading,...
            nodeForce,loadForce,contactForce]...
            =Solver_LoadingMGDCM(obj,mgdcm);
        
        % Solver for DC loading
        [U,UhisLoading,loadHis,strainEnergyLoading,...
            nodeForce,loadForce,contactForce]=Solver_LoadingDC(obj,dc);
        
        % Solver for self folding
        [U,UhisAssemble,strainEnergyAssemble,sprTargetZeroStrain,...
            rotTargetZeroStrain]=Solver_Assemble(obj,selfFold);

        % modify the stiffness matrix for support
        [Kwsupp,T]=Solver_ModKforSupp(obj,K,supp,Tinput,elasticSupportOpen,suppElastic,U);
        
    end    
end