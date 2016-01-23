function moranI = moranI_mesh(vertices_coords,event_vertinds,weight_threshdist)
% moranI_mesh calculates the value of the Moran's I function 
% for events recorded on a mesh, by using the Euclidean distance.
%
% Implemented according to http://cran.r-project.org/web/packages/ape/vignettes/MoranI.pdf
% but with the internal sum over j going only from i+1 till n.
%
% Defined as 
% 	I = \frac{N}{\sum_i \sum_j w_{i,j}} \frac{\sum_i \sum_j w_{i,j}(X_i-X̄)(X_j-X̄)}{\sum_i (X_i - X̄)^2},
% where N is the number of spatial units indexed by i and j (in our case 
% a unit is a vertex of the mesh); X is the variable of interest 
% (whether a vertex was annotated as having an event); X̄ is the mean of X 
% (across all vertices); and wij are spatial weights. A spatial weight w_{ij} 
% was defined to be equal to 1 if the Euclidean distance between i’th and j’th 
% vertices is less than a specific threshold; otherwise, it is equal to 0.
%
% Note that normally the Moran’s I is calculated on a X-Y plane 
% where spatial units are assumed to have the same area and a signal 
% is defined as a quantity of events in a spatial unit. Since this 
% set up is not directly applicable to a 3D mesh, we rather considered
% a spatial unit to be a vertex of a mesh and made the weights w_{ij} 
% to encode the neighborhood property. 
%
% PARAMETERS
% 	vertices_coords: an Nx3-array containing x,y,z coordinates of N vertices of the mesh
% 	event_vertinds: a Nx1 boolean array indicating at which vertices events happened
%	weight_threshdist: a real-valued array of the distance thresholds (can be just one threshold)
%
% Examples for images:
%  1)
%     for k=20:10:100; [XI,YI]=meshgrid(1:k,1:k);Vx=XI(:);Vy=YI(:);Em=zeros(k,k);for i=1:k;for j=1:k;fi=(floor(i/2)*2==i);fj=(floor(j/2)*2==j);Em(i,j)=(fi==fj);end;end;tic;v=moranI_mesh([Vx Vy], Em(:), 1);fprintf('%.0f : %.2f\t',k,v);toc;end
%  2)
%     for k=20:10:100;[XI,YI]=meshgrid(1:k,1:k);Vx=XI(:);Vy=YI(:);Em=[repmat(ones(k,1),1,k/2) repmat(zeros(k,1),1,k/2)];tic;v=moranI_mesh([Vx Vy], Em(:), 1);fprintf('%.0f : %.2f\t',k,v);toc;end
%
% We recommend to apply it to small meshes only (10000 or fewer faces).
% 
% Theodore Alexandrov, EMBL (github theodev)
% Ekaterina Ovchinnikova, KIT (github eovchinn)
% 2015

Nv=size(vertices_coords,1); % number of vertices

xmean=sum(event_vertinds)/Nv; % mean over all vertices with "events"; calc'd this way, since event intensities are asssumed to be 1

summand11=(1-xmean)*(1-xmean); %  when both vertices have events
summand10=(1-xmean)*(0-xmean); % when one vertex has an event
summand00=(0-xmean)*(0-xmean); % when no vertices have events

nomsum=0; % nominator sum
denomsum=0; % nominator sum
Szero=0;

for i=1:Nv % go over all vertices
   Wdist=pdist2( vertices_coords(i+1:end,:), vertices_coords(i,:)  );
   Wi= Wdist <= weight_threshdist;
   Szero=Szero+sum(Wi);
   
   Nj=length(Wi);
   Xi=repmat( event_vertinds(i)-xmean, Nj, 1);
   
   Xj=event_vertinds(i+1:end)-xmean;
   
   nomsum=nomsum + sum(Wi .* Xi .* Xj);
end

denomsum = sum(event_vertinds*summand11 + (~event_vertinds)*summand00);

moranI=Nv/Szero * nomsum/denomsum;

end