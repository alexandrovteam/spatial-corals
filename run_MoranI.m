clear all

%
% This script performs the Moran's I spatial statistical analysis on a 3D model 
% of a coral colony annotated with the GA lesions. 
% 
% A model should be provided as an .stl file. Since Moran's I analysis is super slow,
% we subsampled all models to have 10000 faces (in MeshLab), so the name of the colony
% should contain "_10000faces" e.g. GA_Colony_clipped_4279_10000faces.stl. 
% 
% The annotations should be provided as a .xlsx file, e.g. GA_coordinates_colony_4279.xlsx,
% containing four columns (name of the spot, x, y, z). The first row is the column titles.
% The rows corresponding to lesion spots should have the name starting with 'point'. 
%
% This script requires our implementation of Moran's I function for a 3D mesh.
%
% Theodore Alexandrov, EMBL (github theodev)
% Ekaterina Ovchinnikova, KIT/EMBL (github eovchinn)
% 2014-2016
% 

%------------------ PARAMETERS ------------------

% the folder with stlread.m
addpath('thirdparty'); 

% the folder with 3D models and coordinates of GA lesions
models_folder='/Users/theodore/sci/projects/3D_coral-John/models';

% number of random configurations of spots to be simulated
Nrepeat=100;

% specify the colony ID or several colony IDs to perform analysis on
colonyID='clipped_4279'

% distance thresholds for Moran's I
mI_thresholds=[.01 0.03 0.05 0.07 .1 .2];

%-------------------------------------------------

%
%
% PREPARE THE 3D MODEL
%
% 

tic; % Moran's I is super slow, so we are tracking the time

% load the 3D model
model = stlread([models_folder '/GA_Colony_' colonyID '_10000faces.stl']);

% where the coordinates of lesions will be stored in CSV format after they are loaded from .xlsx
fname_coordsascii=[models_folder '/GA_' colonyID '--lesions_nvcoords.txt'];

% number of vertices
Nv=size(model.vertices,1); 



fprintf('*** Starting processing colony %s... ***\n', colonyID)

% read coordinates of lesion spots
[xlsNUM,xlsTXT] = xlsread([models_folder '/GA_coordinates_colony_' colonyID '.xlsx']);

spots_labels=xlsTXT(2:end,1);
spots_coords=xlsNUM(1:end,1:3);
Nspots=length(spots_labels); % can be diff from Nlesions if few spots have the same nearest vertex

% find lesions spots among all spots (that also include markers and model corners)
lesion_spots_ind=[];
for i=1:size(xlsTXT,1)
    s=xlsTXT{i,1};
    if length(s)>=6 && strcmpi(s(1:5),'point')
        lesion_spots_ind=[lesion_spots_ind; i-1];
    end
end

Nlesions=length(lesion_spots_ind);
lesions_coords=spots_coords(lesion_spots_ind,:);
lesions_labels=spots_labels(lesion_spots_ind);

% for each lesion, find nearest vertex
lesions_vind=false(Nv,1);
for i=1:Nlesions
    nearestv_ind=dsearchn(model.vertices,lesions_coords(i,:));
    lesions_vind(nearestv_ind)=1;
end

%
%
% PERFORM MORAN'S I ANALYSIS
% 
%

fprintf('Calculating Moran''s I for the colony "%s" of %u vertices...\n', colonyID, Nv);

% calculate Moran's I values for the lesion spots for the specific distance thresholds
mIs=zeros(length(mI_thresholds),1);
for it=1:length(mI_thresholds)
    mIs(it)=moranI_mesh(model.vertices, lesions_vind, mI_thresholds(it));
end

%
% simulate random configurations and perform Moran's I analysis for each of them
%

%% simulation of Nrepeat random configurations
rand_mIs=zeros(Nrepeat,length(mI_thresholds));
tmp_allvertind=1:Nv;

% repeat the analysis for Nrepeat random configurations
for ir=1:Nrepeat
    fprlen=fprintf('Simulation %u/%u...', ir, Nrepeat);
    
    % randomly select a subset of vertices from the model whose x-y
    % coordinates are inside the x-y convex hull
    rand_vind=randsample(Nv,Nlesions);
    
    % calc Moran's I for the random configuration for the specified distance thresholds
    for it=1:length(mI_thresholds)
        rand_mIs(ir,it)=moranI_mesh(model.vertices, ismembc(tmp_allvertind,rand_vind)', mI_thresholds(it));
    end
    
    fprintf('%s', repmat(sprintf('\b'),1,fprlen))
    
    % DBG
    %disp('saving all...')
    %save moransI_calcs.mat mI_thresholds mIs rand_mIs
end

%% output values of Moran's I for the lesions and simulated configurations as well as their difference
for it=1:length(mI_thresholds)
    maxMIval_allrandom=max(rand_mIs(:,it));
    fprintf('thresh dist %.2f : real,maxval of simulated,diff(real-maxval) = %g,%g,%g\n', mI_thresholds(it), mIs(it), maxMIval_allrandom,mIs(it)-maxMIval_allrandom)
end


toc % finished tracking time