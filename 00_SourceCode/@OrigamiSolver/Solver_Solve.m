% This function runs the analyses. 
% for some of the analyses, future updates are needed so that code can run.

function Solver_Solve(obj)

    % if the solver is not new and is continuing to solve the status, then
    % automatically initialize the status;
    if obj.continuingLoading==0
        % Initialze the origami structure first
        
        % We assume that the initial deformation field is zero and that the
        % origami has no load applied.
        newNodeNum = size(obj.newNode);
        newNodeNum = newNodeNum(1);
        obj.currentU = zeros(newNodeNum,3);    
        obj.currentAppliedForce = zeros(newNodeNum,3);
        
        % Then calculate the origami mechanical proeprties
        % We first calcualte the folding angle, assume that the current
        % state is stress free. 
        
        if obj.flag2D3D==2
            obj.currentRotZeroStrain=pi*ones(obj.oldCreaseNum,1);
            obj.currentSprZeroStrain=obj.Spr_Theta(...
                obj.currentU,obj.sprIJKL,obj.newNode);
        else  
            if obj.compliantCreaseOpen==1

                obj.compliantCreaseOpen=0;
                obj.Mesh_Mesh();

                obj.currentRotZeroStrain=pi*ones(obj.oldCreaseNum,1);
                obj.currentSprZeroStrain=obj.Spr_Theta(...
                    obj.currentU,obj.sprIJKL,obj.newNode);

                for i = 1:obj.oldCreaseNum
                    if obj.currentSprZeroStrain(i,1) ~=0
                        obj.currentRotZeroStrain(i)=pi-(obj.currentSprZeroStrain(i,1)-pi);
                    end
                end
            % Here we mesh the system with no compliant crease, calculate the
            % rotaiton angle for spring elements and then record the numbers as
            % the rotation angle for Rot.
                obj.compliantCreaseOpen=1;
                obj.Mesh_Mesh();    

            else
                obj.currentRotZeroStrain=pi*ones(obj.oldCreaseNum,1);
                obj.currentSprZeroStrain=obj.Spr_Theta(...
                    obj.currentU,obj.sprIJKL,obj.newNode);

                for i = 1:obj.oldCreaseNum
                    if obj.currentSprZeroStrain(i,1) ~=0
                        obj.currentRotZeroStrain(i)=(obj.currentSprZeroStrain(i,1));
                    end
                end
            end
        
        end
        
        obj.Mesh_MechanicalProperty()
        obj.currentSprZeroStrain=obj.Spr_Theta(...
                obj.currentU,obj.sprIJKL,obj.newNode);

        % We assume that there is no heating applied on to the structure so the
        % system has zero heating input and every node at RT
        obj.currentT=obj.RT*ones(newNodeNum,1);
        obj.currentQ=zeros(obj.envLayer*newNodeNum,1);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Just a note here. The current initialization code is not so
    % perfect. The RotVector and the SprVector sometime is not
    % correctly assocated with each other in some situations. 
    % This needs to be fixed in the future.         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % We sequentially reads the items in the loading controller and
    % perform the analysis. The loading history is recorded. 
    totalLoadSsequance = size(obj.loadingController);
    totalLoadSsequance = totalLoadSsequance(1,2);
    
    for loadSequeceNum=1:totalLoadSsequance
        
        % get the temp loading Controller 
        tempController=obj.loadingController{loadSequeceNum};
        analyzeType=tempController{1};
        tempController=tempController{2};  
        
        if analyzeType=="DC"            
            % perform the displacement controlled loading
            [U,UhisLoading,loadHis,strainEnergyLoading,...
                nodeForce,loadForce,contactForce]=...
                obj.Solver_LoadingDC(tempController); 
            % plot the reults when ploting option is open
            if tempController.plotOpen==1
                obj.Plot_DeformedShape(obj.currentU+obj.newNode, U+obj.newNode);                
            end
            if tempController.videoOpen==1
                obj.Plot_DeformedHis(obj.newNode,UhisLoading); 
            end
            if tempController.detailFigOpen==1                
                obj.Plot_LoadHis(loadHis,UhisLoading)
                obj.Plot_Energy(UhisLoading,strainEnergyLoading);
                obj.Plot_LoadAndReaction(obj.currentU+obj.newNode, ...
                    U+obj.newNode, tempController.load, ...
                    tempController.supp, nodeForce, loadForce);                
                if obj.contactOpen==1
                    obj.Plot_ContactForce(obj.currentU+obj.newNode,...
                    U+obj.newNode,contactForce);
                end
            end             
            % update the current status or origami after loading
            obj.currentAppliedForce=loadForce+obj.currentAppliedForce;  
            obj.currentU=U;
            
            tempController.Uhis=UhisLoading;
            tempController.loadHis=loadHis;
            tempController.strainEnergyHis=strainEnergyLoading;
            
        elseif analyzeType=="NR"
            %perform the displacement controlled loading
            [U,UhisLoading,loadHis,strainEnergyLoading,...
                nodeForce,loadForce,contactForce]=...
                obj.Solver_LoadingNR(tempController);
                        % plot the reults when ploting option is open
            if tempController.plotOpen==1
                obj.Plot_DeformedShape(obj.currentU+obj.newNode, U+obj.newNode);                
            end
            if tempController.videoOpen==1
                obj.Plot_DeformedHis(obj.newNode,UhisLoading); 
            end
            if tempController.detailFigOpen==1                
                obj.Plot_LoadHis(loadHis,UhisLoading)
                obj.Plot_Energy(UhisLoading,strainEnergyLoading);
                obj.Plot_LoadAndReaction(obj.currentU+obj.newNode, ...
                    U+obj.newNode, tempController.load, ...
                    tempController.supp, nodeForce, loadForce);                
                if obj.contactOpen==1
                    obj.Plot_ContactForce(obj.currentU+obj.newNode,...
                    U+obj.newNode,contactForce);
                end
            end             
            % update the current status or origami after loading
            obj.currentAppliedForce=loadForce+obj.currentAppliedForce;  
            obj.currentU=U;
            
            tempController.Uhis=UhisLoading;
            tempController.loadHis=loadHis;
            tempController.strainEnergyHis=strainEnergyLoading;
            
        elseif analyzeType=="MGDCM"
            % perform the MGDCM controlled loading
            [U,UhisLoading,loadHis,strainEnergyLoading,...
                nodeForce,loadForce,contactForce]...
                =Solver_LoadingMGDCM(obj,tempController);
            if tempController.plotOpen==1
                obj.Plot_DeformedShape(obj.currentU+obj.newNode, U+obj.newNode);                
            end
            if tempController.videoOpen==1
                obj.Plot_DeformedHis(obj.newNode,UhisLoading); 
            end
            if tempController.detailFigOpen==1                
                obj.Plot_LoadHis(loadHis,UhisLoading)
                obj.Plot_Energy(UhisLoading,strainEnergyLoading);
                obj.Plot_LoadAndReaction(obj.currentU+obj.newNode, ...
                    U+obj.newNode, tempController.load, ...
                    tempController.supp, nodeForce, loadForce);                
                if obj.contactOpen==1
                    obj.Plot_ContactForce(obj.currentU+obj.newNode,...
                    U+obj.newNode,contactForce);
                end
            end            
            % update the current status or origami after loading
            obj.currentAppliedForce=loadForce+obj.currentAppliedForce;  
            obj.currentU=U;
            
            tempController.Uhis=UhisLoading;
            tempController.loadHis=loadHis;
            tempController.strainEnergyHis=strainEnergyLoading;
            
        elseif analyzeType=="SelfFold"             
            % perform the self folding analysis
            [U,UhisAssemble,strainEnergyAssemble,...
                sprTargetZeroStrain,rotTargetZeroStrain]=...
                obj.Solver_Assemble(tempController);
            % plot the reults when ploting option is open
            if tempController.plotOpen==1
                obj.Plot_DeformedShape(obj.currentU+obj.newNode, U+obj.newNode);
            end
            if tempController.videoOpen==1
                obj.Plot_DeformedHis(obj.newNode,UhisAssemble); 
            end
            if tempController.detailFigOpen==1                
                obj.Plot_Energy(UhisAssemble,strainEnergyAssemble);
            end
            % update the current status or origami after loading
            obj.currentU=U;
            obj.currentRotZeroStrain=rotTargetZeroStrain;
            obj.currentSprZeroStrain=sprTargetZeroStrain;
            
            tempController.Uhis=UhisAssemble;
            tempController.strainEnergyHis=strainEnergyAssemble;
            
        elseif analyzeType=="ThermalLoading"
            % perform the thermal loading analysis
            [U,UhisThermal,energyHisThermal,temperatureHistory,...
                rotTargetZeroStrain,sprTargetZeroStrain]=...
                obj.Solver_LoadingThermal(tempController);
            
            % plot the reults when ploting option is open
            if tempController.plotOpen==1
                obj.Plot_DeformedShapeTemp(tempController,...
                    obj.currentU+obj.newNode, U+obj.newNode, ...
                    squeeze(temperatureHistory(:,tempController.thermalStep)));
            end
            if tempController.videoOpen==1
                obj.Plot_DeformedHisTemp(obj.newNode,UhisThermal,temperatureHistory); 
            end
            if tempController.detailFigOpen==1                
                obj.Plot_Energy(UhisAssemble,energyHisThermal);
            end
            % update the current status or origami after loading
            obj.currentU=U;
            obj.currentT=temperatureHistory(:,tempController.thermalStep);
            obj.currentRotZeroStrain=rotTargetZeroStrain;
            obj.currentSprZeroStrain=sprTargetZeroStrain;
            
            tempController.Uhis=UhisThermal;
            tempController.strainEnergyHis=energyHisThermal;
            tempController.temperatureHis=temperatureHistory;
            
        end
    end
end