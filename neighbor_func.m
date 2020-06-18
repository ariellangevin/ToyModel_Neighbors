function [dxdt, x, map] = neighbor_func(t,x,abDose,map,Kout,abType) %CONFIG)

% INPUTS
% abDose = antibiotic concentration
% map = number of cells / configuration
% Kout = type of cell in map
% abType = type of antibiotic

% OUTPUTS
% x(odd) = concentrations levels for each pos in map
% x(even) = biomass of cells for each pos in map
% map = spatial arrangement of cells

% Example call:
% [t,x] = ode45(@neighbor_func,[0 120],ones(1,73),[],'WT','DacrB',1,2);

% Originally written by A.M. Langevin, last edited 6/6/2018 by A.M. Langevin

%% Initialization

area=size(map);
dim=area(1);
Kin=1;

mapId=reshape(1:dim^2,dim,dim);

% we can do the same for growth rate too

linMap=reshape(map,[],1);
linKout=reshape(Kout,[],1); % reshape(A,1,[]) % reshape matrix into array for Kin and Kout so these can
% easily be passed along when cell divides

% Initialize solutions keeping in mind the total number of possible cells

% Step into CA here to arrange

for i=1:length(linMap)
    N_map_Id=[zeros(1,dim);mapId(1:(dim-1),:)];
    S_map_Id=[mapId(2:dim,:);zeros(1,dim)];
    E_map_Id=[mapId(:,2:dim),zeros(dim,1)];        
    W_map_Id=[zeros(dim,1),mapId(:,1:(dim-1))];

    Coor_map_Id(i,:)=[N_map_Id(i) S_map_Id(i) E_map_Id(i) W_map_Id(i)];
    
    N_map=[zeros(1,dim);map(1:(dim-1),:)];
    S_map=[map(2:dim,:);zeros(1,dim)];
    E_map=[map(:,2:dim),zeros(dim,1)];        
    W_map=[zeros(dim,1),map(:,1:(dim-1))];

    Coor_map(i,:)=[N_map(i) S_map(i) E_map(i) W_map(i)];
    
    % find which ones have cells and leave the proper address, but update
    % those who have no cells in assigned position to zero
    
    mapNeighbors=Coor_map.*Coor_map_Id;
    
end

% Neighbor Associations

neighbors=2*mapNeighbors+1;
neighbors(neighbors==1)=2*length(linMap)+1;

cellsAt=2*find(linMap==1)-1;
%if abDose==0
    %x(2*length(linMap)+1)=abDose;
%else
    x(2*length(linMap)+1)=abDose-sum(x(cellsAt))/(dim^2-3);
%end

KoutAll=Kout;
KoutAll(KoutAll==0)=Kin;
    
% identify who is neighboring who for the diffEqs
% make a list of neighbors

%% Differential equations

% Specify constants

%mu = 0.0116/1.5; % grows half as fast on MGC media % in lb mu=0.19 for Dacrb mu=0.16 for WT with pumps

if strcmpi(abType,'Cm')==1
    Kc = 0.143; % MIC of cells without pumps
    hc = 1.43; 
    mu = 0.0116; %doubling time is 1 h
elseif strcmpi(abType,'Cp')==1
    Kc = 0.890; % MIC of cells without pumps
    hc = 1.80;
    mu = 0.0116/1.5; %doubling time is 1.5h
end

% Run differential equations

%x(2*length(linMap)+1)=abDose;

for i=1:2:2*length(linMap)

    j=round(i/2);
    
    % odd numbers correspond to intracellular concentration of cells
    % even numbers correspond to the biomass of cells
    % the last number corresponds to the extracellular concentration
    
    if linMap(j)==0
        x(i)=x(2*length(linMap)+1);
        dxdt(i)=0;
        dxdt(i+1)=0;
    end
    
    if linMap(j)==1
%         if linKout(j)>1 % from Xi's data
%             mu=0.0773/10;
%         elseif linKout(j)==1
%             mu=0.1107/10;
%         end
        dxdt(i)=2*((sum(1/8*((KoutAll(Coor_map_Id(j,:)))+1).*x(neighbors(j,:))'./2)+0.5*x(2*length(linMap)+1)-Kout(j).*x(i)));
        dxdt(i+1)=mu*x(i+1)*1/(1 + (x(i)/Kc)^hc);
    end
    
    % Optional code for dividing cells
    
        % if there is a cell there, check if positions need updating,
        % NO DIVISION, if main cell reaches given biomass, stop experiment
        % then and calculate growth from time experiment ran instead
     
    
%     CheckHood=find(mapNeighbors(j,:)==0);
    
%     if x(i+1)>=2 && isempty(CheckHood)==0
%         % make the new cell, divide the biomass and Ab dose
%         %CheckHood=find(mapNeighbors(j,:)==0);
%         r=randi(length(CheckHood));
%         NewCell=Coor_map_Id(j,r);
%         map(NewCell)=1;
%         Kout(NewCell)=Kout(j);
%         x(NewCell*2+1)=x(i)/2;
%         x(NewCell*2+2)=1;
%         x(i+1)=x(i+1)-1;
%         x(i)=x(i)/2;
%         linMap=reshape(map,[],1); % update linMap
% %     elseif x(i+1)>=2 && isempty(CheckHood)==1
% %         % move neighbor
% %         r=randi(4);
% %         NeighborToMove=mapNeighbors(j,r);
% %         CheckNewHood=find(mapNeighbors(NeighborToMove,:)==0);
% %         r=randi(length(CheckNewHood));
% %         MoveTo=Coor_map_Id(NeighborToMove,r);
% %         map(MoveTo)=1;
% %         Kout(MoveTo)=Kout(NeighborToMove);
% %         Kout(NeighborToMove)=Kout(j);
% %         x(MoveTo*2+1)=x(NeighborToMove*2+1);
% %         x(MoveTo*2+2)=x(NeighborToMove*2+2);
% %         x(NeighborToMove*2+1)=x(i)/2;
% %         x(NeighborToMove*2+2)=1;
% %         x(i+1)=x(i+1)-1;
% %         x(i)=x(i)/2;
% %         linMap=reshape(map,[],1); % update linMap
%     end
    
end

dxdt(2*length(linMap)+1)=0; % don't accumulate biomass if no cell is there
% extracellular concentration of antibiotic, assume instant diffusion

% stop at certain time points for CA snap shot and if the cells divid

dxdt = dxdt';
