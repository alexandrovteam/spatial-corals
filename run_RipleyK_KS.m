clear all

%------------------ PARAMETERS ------------------

% the folder with stlread.m
addpath('thirdparty'); 

% the folder with 3D models and coordinates of GA lesions
models_folder='/Users/theodore/sci/projects/3D_coral-John/models';

% number of random configurations of spots to be simulated
Nrepeat=1000;

% specify the colony ID or several colony IDs to perform analysis on
colonyID='clipped_4279'

% whether to use the edge correction for Ripley's K
RipleyK_edgecorr=0;

%-------------------------------------------------

%
%
% PREPARE THE 3D MODEL AND RENDER IT
%
% 

% load the 3D model
model = stlread([models_folder '/GA_Colony_' colonyID '.stl']);

% 
fname_coordsascii=[models_folder '/GA_' colonyID '--lesions_nvcoords.txt'];
% fname_coordsascii='/Users/theodore/sci/projects/3D_coral-John/models/GA_6085--lesions_nvcoords.txt';

% number of vertices
Nv=size(model.vertices,1); 


%%
fprintf('*** Starting processing colony %s... ***\n', colonyID)

%% read coordinates of lesion spots
[xlsNUM,xlsTXT] = xlsread([models_folder '/GA_coordinates_colony_' colonyID '.xlsx']);

spots_labels=xlsTXT(2:end,1);
spots_coords=xlsNUM(1:end,1:3);
Nspots=length(spots_labels);

%% find lesions spots among all spots (that also include markers and model corners)
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

%% for each lesion, find nearest vertex
lesions_nvcoords=zeros(Nlesions,3);
for i=1:Nlesions
    nearestv_ind=dsearchn(model.vertices,lesions_coords(i,:));
    lesions_nvcoords(i,:)=model.vertices(nearestv_ind,:);
end

%% save the coords of the lesions into an ascii file
fid=fopen(fname_coordsascii,'w');
for i=1:size(lesions_nvcoords)
    fprintf(fid, '%.2f,%.2f,%.2f\n', lesions_nvcoords(i,1), lesions_nvcoords(i,2), lesions_nvcoords(i,3));
end
fclose(fid);

%% prepare a subset of vertices which are inside the convex hull x-y projection of the lesion spots
%%% *** not needed anymore, because John clipped the colonies ***
%%%chull_sind=convhull(lesions_coords(:,1),lesions_coords(:,2)); 
%%%hull_x=[lesions_coords(chull_sind,1);lesions_coords(chull_sind(1),1)]; % x coords of the hull vertices
%%%hull_y=[lesions_coords(chull_sind,2);lesions_coords(chull_sind(1),2)]; % y coords of the hull vertices
%%%vertices_inhull=inpolygon(model.vertices(:,1),model.vertices(:,2),hull_x,hull_y);
%%%inhull_vind=find(vertices_inhull);

% remnants of using chull for not-clipped colonies (ideally should be refactored)
vertices_inhull=true(size(model.vertices,1),1);
inhull_vind=(1:size(model.vertices,1))';
Nv_insidehull=length(inhull_vind); % number of vertices inside the hull

%% visualisation of x-y of the lesion spots and the convex hull
figure;
scatter(lesions_coords(:,1),lesions_coords(:,2)); 
title('X-y coordinates of the lesion spots')
xlabel('x');ylabel('y');
saveas(gcf, [colonyID '_lesions.png'])



% create a new model structure, only inside the lesions area
model_lesionsarea=struct('faces',[],'vertices',model.vertices);
for i=1:size(model.faces,1)
    %fprlen=fprintf('%u',i);
    
    vi1=model.faces(i,1);
    vi2=model.faces(i,2);
    vi3=model.faces(i,3);
    
    if vertices_inhull(vi1) && vertices_inhull(vi2) && vertices_inhull(vi3)
        model_lesionsarea.faces=[model_lesionsarea.faces; model.faces(i,:)];
    end
    
    %fprintf('%s', repmat(sprintf('\b'),1,fprlen))
end

% render the 3D model
figure;clf
set(gca,'FontName','Arial','FontSize',14)
patch(model_lesionsarea,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);
xlabel('x','FontName','Arial','FontSize',16)
ylabel('y','FontName','Arial','FontSize',16)
zlabel('z','FontName','Arial','FontSize',16)
grid

view([0 0 1])

% indicate lesion spots with balls
radius=0.02;
xoff=0.01;
yoff=0.01;
zoff=0.03;
%colors=colormap(jet(Nlesions));
for i=1:Nlesions
    hold on;
    
    xc=lesions_coords(i,1);
    yc=lesions_coords(i,2);
    zc=lesions_coords(i,3);
    
    [ellx,elly,ellz]=ellipsoid(xc,yc,zc,radius,radius,radius,10);
    %he=surf(ellx,elly,ellz,'LineStyle','none','FaceColor',colors(i,:));
    he=surf(ellx,elly,ellz,'LineStyle','none','FaceColor',[255 153 0]/255);
    %text(xc+xoff,yc+yoff,zc+zoff,sprintf('%s',spots_labels{i}),'FontSize',8,'Interpreter','non');
end

str=sprintf('Colony "%s", the area including GA lesions', strrep(colonyID,'clipped_',''));
title(str,'FontName','Arial','FontSize',18)

% save the view of the 3D model with lesion spots indicated as a PNG
saveas(gcf, [colonyID '-3D_model_lesions.png'])

%
%
% PERFORM RIPLEY'S K AND KOLMOGOROV-SMIRNOV ANALYSES
% 
%

% pairwise distances between the lesion spots
real_dists=pdist(lesions_nvcoords);

% prepare the boundbox for the Ripley's K analysis
hmin=0;
hmax=max([ max(real_dists) max(model.vertices(:,1))-min(model.vertices(:,1)) max(model.vertices(:,2))-min(model.vertices(:,2)) ]);
N_hsteps=20;
hstep=(hmax-hmin)/(N_hsteps-1);
hvals=hmin:hstep:hmax;

boundbox=[min(lesions_nvcoords(:,1)) max(lesions_nvcoords(:,1)) min(lesions_nvcoords(:,2)) max(lesions_nvcoords(:,2))];

% calculate Ripley's K function for the lesion spots
real_rK=ripleyK(lesions_nvcoords(:,1:2),hvals,boundbox, RipleyK_edgecorr);

%
% simulate random configurations and perform Ripley's K and Kolmogorov-Smirnov analyses for each of them
%

% storage for the statistics for the simulated configurations
pvals=zeros(Nrepeat,1);
ksvals=zeros(Nrepeat,1);
rK_vals=zeros(Nrepeat,N_hsteps);

wb = waitbar(0,'Simulating random configurations & calculating statistis on them...');

%repeat random sampling
for ir=1:Nrepeat 
     waitbar(ir/Nrepeat,wb);
     
    fprlen=fprintf('%u',ir);
    
    % randomly select a subset of vertices from the model whose x-y
    % coordinates are inside the x-y convex hull
    rand_vind=randsample(Nv_insidehull,Nlesions);

    % get the coordinates of randomly sampled spots
    rand_vcoords=model_lesionsarea.vertices(inhull_vind(rand_vind),:);
    
    % distances between randomly sampled spots
    rand_dists=pdist(rand_vcoords);
    
    % Kolm-Smirn test to test for difference
    [~,pvals(ir),ksvals(ir)]=kstest2(real_dists,rand_dists);
    
    % Ripley K function for the randomly sampled spots
    rK_vals(ir,:)=ripleyK(rand_vcoords(:,1:2),hvals,boundbox, RipleyK_edgecorr);
    
    fprintf('%s', repmat(sprintf('\b'),1,fprlen))
end
close (wb);

% calculate values of Ripley's K for tightly clustered data
clust_vind=1:floor(Nv_insidehull/10):Nv_insidehull-100;
while length(clust_vind)<Nlesions
    clust_vind=[clust_vind clust_vind+1];
end
clust_vind=sort( clust_vind(1:Nlesions) );
clust_vcoords=model_lesionsarea.vertices(inhull_vind(clust_vind),:);
rK_clust=ripleyK(clust_vcoords(:,1:2),hvals,boundbox, RipleyK_edgecorr);

% calculate values of Ripley's K for uniformly spaced data
uniform_vind=1: floor(Nv/Nlesions) :Nv;
uniform_vind=uniform_vind(1:Nlesions);
uniform_vcoords=model_lesionsarea.vertices(inhull_vind(uniform_vind),:);
rK_uniform=ripleyK(uniform_vcoords(:,1:2),hvals,boundbox, RipleyK_edgecorr);

%% statistics of difference between Ripley's K of random configurations and the lesion spots
rK_quantlow=quantile(rK_vals,  0.0,  1);
rK_quanthigh=quantile(rK_vals,  1.0,  1);

diff_real_quanthigh=max(real_rK-rK_quanthigh, 0);
diff_real_quantlow=max(rK_quantlow-real_rK, 0);

[valH,indH]=max(diff_real_quanthigh);
[valL,indL]=max(diff_real_quantlow);

if valH>valL
    rK_maxdiffH=true; % achieved on highquantile?
    rK_maxdiff=valH/max(real_rK); % value
    rK_maxdiffind=indH; % index where achieved
else
    rK_maxdiffH=false;
    rK_maxdiff=valL/max(real_rK);
    rK_maxdiffind=indL;
end
    
%
%
% RENDER AND SAVE PLOTS FOR RIPLEY'S K AND KOLMOGOROV-SMIRNOV ANALYSES
% 
%
 

% plots for Ripley's K
figure;
plot(hvals,real_rK,'color','r')
hold on; plot(hvals,rK_vals(1,:),'color',[.5 .5 .5]);
hold on;plot(hvals,rK_quantlow,'color','b', 'LineWidth',2)
hold on;plot(hvals,rK_quanthigh,'color','g', 'LineWidth',2)
for i=1:Nrepeat
    plot(hvals,rK_vals(i,:),'color',[.5 .5 .5]);
    hold on;
end
hold on;plot(hvals,rK_quantlow, 'color','g', 'LineWidth',2)
hold on;plot(hvals,rK_quanthigh, 'color','b', 'LineWidth',2)
hold on;plot(hvals,real_rK, 'color','r', 'LineWidth',2)

le=legend('GA lesions','all simulated configurations','maximal simulated, Kmin','minimal simulated, Kmax', 'Location','southeast');
set(le,'FontName','Arial','FontSize',16)

if rK_maxdiffH
    line([hvals(rK_maxdiffind) hvals(rK_maxdiffind)],[real_rK(rK_maxdiffind) rK_quanthigh(rK_maxdiffind)],'Color','m')
    text(hvals(rK_maxdiffind)*1.05, rK_quanthigh(rK_maxdiffind)+abs(real_rK(rK_maxdiffind)-rK_quanthigh(rK_maxdiffind))/2,'maxdiff','Color','m','FontSize',16,'FontName','Arial')
else
    line([hvals(rK_maxdiffind) hvals(rK_maxdiffind)],[real_rK(rK_maxdiffind) rK_quantlow(rK_maxdiffind)],'Color','m')
    text(hvals(rK_maxdiffind)*1.05, rK_quantlow(rK_maxdiffind)-abs(real_rK(rK_maxdiffind)-rK_quantlow(rK_maxdiffind))/2,'maxdiffK','Color','m','FontSize',16,'FontName','Arial')
end

ylim([min(real_rK) max(real_rK)*1.05])
set(gca,'FontName','Arial','FontSize',16)
xlabel('search radius, t', 'FontSize', 20)
ylabel('Ripley''s K', 'FontSize',20)

saveas(gcf, [colonyID '-RipleyK_statistics.png'])


% plots for the KS test
figure;
plot(pvals)
title({'p-values of the Kolmogorov-Smirnov test', ...
sprintf('mean p-value=%.1g, median p-value=%.1g, %.1f%% of p-values are <0.01/%u\n', ...
               mean(pvals), median(pvals), sum(pvals<0.01/Nrepeat)/Nrepeat*100, Nrepeat)})
xlabel('Random configurations of spots')

saveas(gcf, [colonyID '-KS_test_pvalues.png'])

figure;
subplot(211);
hist(real_dists);
title('Histogram of distances between lesions');
xlabel('Distance between two lesions');
ylabel('Frequency')

subplot(212);
hist(rand_dists);
title('Histogram of distances between randomly distributed spots (exemplary random configuration of spots)');
xlabel('Distance between two random spots');
ylabel('Frequency')

saveas(gcf, [colonyID '-histograms_distances.png'])


% output the results in console
fprintf('Colony name,mean p-value,median p-value,ratio of p-values <0.01/%u,max diff\n', Nrepeat);
fprintf('%s,%.1g,%.1g,%.0f%%,%.2f\n',colonyID, mean(pvals), median(pvals), sum(pvals<0.01/Nrepeat)/Nrepeat*100,rK_maxdiff);



%% statistics for baseline
% real_dists=rand_dists;
% 
% pvals=zeros(Nrepeat,1);
% ksvals=zeros(Nrepeat,1);
% for ir=1:Nrepeat %repeat random sampling
%     fprlen=fprintf('%u',ir);
%     
%     rand_vind=randsample(Nv_insidehull,Nlesions);
%     rand_vcoords=model_lesionsarea.vertices(inhull_vind(rand_vind),:);
%     rand_dists=pdist(rand_vcoords);
%     
%     [~,pvals(ir),ksvals(ir)]=kstest2(real_dists,rand_dists);
%     
%     fprintf('%s', repmat(sprintf('\b'),1,fprlen))
% end
% 
% %%
% figure;
% plot(pvals)
% title({'p-values of the Kolmogorov-Smirnov test (random-random)', ...
%        sprintf('mean p-value=%.1g, median p-value=%.1g, %.1f%% of p-values are <0.01/%u\n', ...
%                mean(pvals), median(pvals), sum(pvals<0.01/Nrepeat)/Nrepeat*100, Nrepeat)})
% xlabel('Random configurations of spots')
% 
% saveas(gcf, [colonyID '-KS_test_pvalues_baseline.png'])
% 
% figure;
% subplot(211);
% hist(real_dists);
% title('Histogram of distances between randomly distributed spots (exemplary random configuration of spots)');
% xlabel('Distance between two lesions');
% ylabel('Frequency')
% 
% subplot(212);
% hist(rand_dists);
% title('Histogram of distances between randomly distributed spots (exemplary random configuration of spots)');
% xlabel('Distance between two random spots');
% ylabel('Frequency')
% 
% saveas(gcf, [colonyID '-histograms_distances_baseline.png'])
% 

