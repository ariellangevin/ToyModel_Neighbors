% Hardcode for Neighbors

% This is a toy model describing how export of antibiotics by certain
% cells can produce harmful effects to their neighbors. These scripts simulate
% cell growth with different neighbors in the presence of antibiotics, using
% an agent-based model with Moore neighborhood architecture to represent the
% spatial interactions between cells and the environment. Each cell's growth
% and antibiotic concentration are represented by ordinary differential
% equations. Please see associated manuscript below for more details.

% Originally written by A.M. Langevin, last edited 6/6/2018 by A.M. Langevin

% Related manuscript: Wen, X., Langevin, A.M. & Dunlop, M.J. Antibiotic export
% by efflux pumps affects growth of neighboring bacteria. Sci Rep 8, 15120
% (2018). https://doi.org/10.1038/s41598-018-33275-4.

close all, clc, clear

% Select Main and Select abType

Configuration=2; % can look at more than one configuration
    % see lines 66-76, optional configurations could be further expanded
Neighbor={'WT','DacrB'}; % which neighbors the main cell has
Main={'WT','DacrB'}; % which cell type to explore
abTypes={'Cm','Cp'}; % which antibiotic is 

n=0; % initialization

for mA=1:length(Main)
    cellTypeMain=Main(mA);
    for mB=1:length(abTypes)
        abType=abTypes(mB);

        % how long to run experiment (dividing time 60 min)
        Tfinal=60;

        % antibiotic concentrations to assess for each antibiotic
            % these are just set from Xi et al., but can vary
        
        if strcmpi(abType,'Cm')==1
            abDosage=[0 0.1 0.2 0.5 1 2];
        elseif strcmpi(abType,'Cp')==1
            abDosage=[0 0.02 0.1 0.2 0.3 0.6];
        end

        % Look at configurations
        
        for i=1:length(Configuration)
            CONFIG=Configuration(i);

            for j=1:length(Neighbor)
                cellTypeNeighbor=Neighbor(j);

                for k=1:length(abDosage)
                    abDose=abDosage(k);
                    %cellTypeMain=Main;

                    % Need to configure the initial map outside or else Matlab gets confused

                    % Initialization of configuration

                    dim=5; % size of map

                    map=zeros(dim); % hard clear map

                    if CONFIG==0
                        startMap=[0 0 0; 0 1 0; 0 0 0];
                    elseif CONFIG==1
                        startMap=[0 0 0; 0 1 1; 0 0 0];
                    elseif CONFIG==2
                        startMap=[0 1 0; 0 1 1; 0 0 0];
                    elseif CONFIG==3
                        startMap=[0 1 0; 0 1 1; 0 1 0];
                    elseif CONFIG==4
                        startMap=[0 1 0; 1 1 1; 0 1 0];
                    end

                    map(round(dim/2)-1:round(dim/2)+1,round(dim/2)-1:round(dim/2)+1)=startMap; % center cells
                        % this maps is to be used for cells to divide - this was
                        % not done in Xi et al. but would be a very cool adaption

                    % Hardcoded constants

                    % Influx and efflux of antibioitcs

                    Kin=1;
                    if strcmpi(abType,'Cm')==1
                        Kout_WT=4;
                    elseif strcmpi(abType,'Cp')==1
                        Kout_WT=3;
                    end
                    Kout_DacrB=Kin;

                    % Influx and efflux of neighboring cells

                    if strcmpi(cellTypeNeighbor,'WT')==1
                        KoutNeighbor=Kout_WT;
                    elseif strcmpi(cellTypeNeighbor,'DacrB')==1
                        KoutNeighbor=Kout_DacrB;
                    end

                    % Influx and efflux of main cell based on neighboring cells

                    if strcmpi(cellTypeMain,cellTypeNeighbor)==1
                        Kout=KoutNeighbor*map;
                    elseif strcmpi(cellTypeMain,'WT')==1 && strcmpi(cellTypeMain,cellTypeNeighbor)==0
                        Kout=KoutNeighbor*map;
                        Kout(round(dim/2),round(dim/2))=Kout_WT;
                    elseif strcmpi(cellTypeMain,'DacrB')==1 && strcmpi(cellTypeMain,cellTypeNeighbor)==0
                        Kout=KoutNeighbor*map;
                        Kout(round(dim/2),round(dim/2))=Kout_DacrB;
                    end

                    % Itialize cells

                    x0=[zeros(1,dim^2);ones(1,dim^2)];
                    x0_up=[reshape(x0,1,[]),0];

                    linMap=reshape(map,[],1);
                    for y=1:2:2*length(linMap)
                        z=round(y/2);
                        if linMap(z)==0
                            x0_up(y)=abDose;
                            x0_up(y+1)=0;
                        end
                    end

                    % Run simulation

                    [t,x] = ode45(@neighbor_func,[0 Tfinal],x0_up,[],abDose,map,Kout,abType);

                    % Find relevant cells and sort by neighbors or not

                    mapId=reshape(1:dim^2,dim,dim);

                    MainLoco=mapId(round(dim/2),round(dim/2));
                    NeighborIdx=find(Kout~=0);
                    NeighborIdxNoMain=NeighborIdx(find(NeighborIdx~=MainLoco));

                    MainCell=2*MainLoco;
                    NeighborCells=2*NeighborIdxNoMain;

                    % Save compiled data

                    All_Main(i,j,k)=(x(length(t),MainCell)-x(1,MainCell))/t(length(t));
                    All_Neighbors(i,j,k)=nanmean((x(length(t),NeighborCells)-x(1,NeighborCells))/t(length(t)));
                    Std_Neighbors(i,j,k)=nanstd((x(length(t),NeighborCells)-x(1,NeighborCells))/t(length(t)));


                    % Optional Figures of individual configurations

        %             linesizing=1;
        %            
        %             figure(n+5)
        %             subplot(2,1,1)
        %             if strcmpi(cellTypeMain,'DacrB')==1
        %                 plot(t,x(:,MainCell),'m-','LineWidth',linesizing)
        %             elseif strcmpi(cellTypeMain,'WT')==1
        %                 plot(t,x(:,MainCell),'g-','LineWidth',linesizing)
        %             end
        %             hold on
        %             for p=1:length(NeighborCells)
        %                 Neigh=NeighborCells(p);
        %             if strcmpi(cellTypeNeighbor,'DacrB')==1
        %                 plot(t,x(:,Neigh),'m-.')
        %             elseif strcmpi(cellTypeNeighbor,'WT')==1
        %                 plot(t,x(:,Neigh),'g-.')
        %             end
        %             end
        %             plot(t,ones(length(t),1),'k--')
        %             plot(t,ones(length(t),1)+1,'k--')
        %             hold off
        %             title([num2str(abDose),' µg/mL Cm'])
        %             %title('Biomass of cells')
        %             legend('Main Cell','Neighboring Cells','location','northwest')
        %             axis([0 Tfinal 0 3])
        %             ylabel('Growth Rate (min^{-1})')
        %             %xlabel('Time, min')
        %             set(gca,'fontsize', 14)
        % 
        %             figure(n+5)
        %             subplot(2,1,2)
        %             if strcmpi(cellTypeMain,'DacrB')==1
        %                 plot(t,x(:,MainCell-1),'m-','LineWidth',linesizing)
        %             elseif strcmpi(cellTypeMain,'WT')==1
        %                 plot(t,x(:,MainCell-1),'g-','LineWidth',linesizing)
        %             end
        %             hold on
        %             for l=1:length(NeighborCells)
        %                 Neigh=NeighborCells(l)-1;
        %             if strcmpi(cellTypeNeighbor,'DacrB')==1
        %                 plot(t,x(:,Neigh),'m-.')
        %             elseif strcmpi(cellTypeNeighbor,'WT')==1
        %                 plot(t,x(:,Neigh),'g-.')
        %             end
        %             end
        %             %plot(t,ones(length(t),1)*abDose,'k--')
        %             hold off
        %             %title('Intracellular concentration of cells')
        %             %legend('Main Cell','Neighboring Cells','location','northwest')
        %             axis([0 Tfinal 0 1.5])
        %             xlabel('Time, min')
        %             ylabel('[Chloramphenicol]_{int}, µg/mL')
        %             set(gca,'fontsize', 14)

                    % update these for when time is varying
                    n=n+1;

                end
            end
        end


        % Compiling data and visualizing experimental results from different configurations

        for i=1:length(Configuration)

            % Compiling results from experiments
            MainResults=[];
            for k2=1:length(abDosage)
                MainResults=[MainResults,All_Main(i,:,k2)];
            end

            %NeighborResults=[All_Neighbors(i,:,1),All_Neighbors(i,:,2)];

            %MainError=[zeros(1,4)]; % only 1 sample, change to Nan if problematic
            %NeighborError=[Std_Neighbors(i,:,1),Std_Neighbors(i,:,2)];

            %Z=[NeighborResults; MainResults];
            %errZ=[NeighborError; MainError];

            % READ THIS for understanding x-axis:
            % Yaxis 1=Neighboring Cells, 2= Main Cells
            % Xaxis 1= Wt with WT neighbors, 2=DacrB cells with WT neighbors, 3=Wt with DacrB neighbors, 4=DacrB cells with DacrB neighbors
            % % 
            figure(mA+2*(mB-1))
            %subplot(length(Main),length(abTypes),(mA+2*(mB-1)))
            p=bar(MainResults,'grouped')
            xlabel(abType)
            ylabel({'Growth Rate (min^{-1})'})
            title([cellTypeMain,' cells'])
            %title({'Effect of neighboring cells on growth rate';['Results for ',num2str(Configuration(i)),' neighbors']})
            xticklabels({'WT','\Delta acrB','WT','\Delta acrB','WT','\Delta acrB','WT','\Delta acrB','WT','\Delta acrB','WT','\Delta acrB'})
            %xticklabels({'WT : WT','\Delta acrB : WT','WT : \Delta acrB','\Delta acrB : \Delta acrB'})
            set(gca,'fontsize', 14)
            set(p,'FaceColor','blue')
            %set(gca,'yScale','log')
            axis([0 13 0 0.025])
            % 

        end

    end
end


