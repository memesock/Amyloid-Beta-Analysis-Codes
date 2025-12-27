
function SphagettiPlot(folderName, species, seq, RMSD_binary, G, clusters, spatialVarUpperLim, PARtrueOrder)
%% Plot spectrally clustered conformers
%   Requires 'classAverageDisordered' to be run first, generating a
%   'clusters' variable. 
%
%   RNA backbone searh uses phosphorus point sampling method from tortuosity scripts  
%
%   GW - 2023 March
%      - updated 2024 January to run without PyMol dependencies
%


%% Create and populate folders with previously identified clusters
%Delete the previous clustering results to avoid confusion in case you generated new clusters
if exist([folderName,'/PDBs_SpectralClustered'],'dir')
    rmdir([folderName,'/PDBs_SpectralClustered'],'s')
end

[files_grouped,nClusters,clusterPDBnumbers] = GroupClusters(folderName, clusters);

%% Define subplot layout 
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(4,floor(nClusters/2)+1)

%% Plot graph from classAverageDisordered.m
nexttile([2 1])
G = graph(RMSD_binary,'omitselfloops','upper');

h = plot(G,'-db','LineWidth',1,'MarkerSize',5);
set(gcf,'color','w')
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
box on
grid off
h.EdgeColor = [0.5 0.5 0.5];
h.Marker = 'O';

% Color graph nodes by cluster
colors = colormap(jet);
colorSpacing = floor(numel(colors(:,1)) / (nClusters+1) * 0.9);
for j = 1:nClusters+1
    colorsSpaced(j,:) = colors(j*colorSpacing,:);

    clusterIndices{j} = clusters(clusters(:,2)==j-1);
    highlight(h,clusterIndices{j},'NodeColor',colorsSpaced(j,:))
end


%% Perform secondary alignment of structures within each spectral cluster 
if ~exist([folderName,'/PDBs_SpectralClustered/Cluster0/Aligned'],'file') %Don't repeat this if already done
    alignClusters_allToOne(folderName,files_grouped,nClusters,clusterPDBnumbers)
end

disp('Generating and saving figures...')
%% Sample backbone coordinates thru species-dependent coarse graining

for clust = 0:(nClusters-1)

    files_inThisCluster = files_grouped{clust+1};
    nStructures = numel(files_inThisCluster);

    figure; hold all
    title(['Cluster ',num2str(clust),'; # of structures = ',num2str(nStructures)],'Color',colorsSpaced(clust+1,:))
    set(gcf,'color','w')

    for name = 1:nStructures

        pdb = pdbread([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/Aligned/',files_inThisCluster{name}]);
        %pdb = pdbread([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/',files_inThisCluster{name}]);
        pdb = pdb.Model.Atom;

        [S] = getWireframeIndividual(pdb,species,PARtrueOrder);

        plot3(S(:,1),S(:,2),S(:,3),'.-','Color',[0.95 0.95 0.95]) % Plot wireframe of each structure with nonzero transparency
       
        % Save coordinates for downstream averaging
        S_allX(:,name) = S(:,1);
        S_allY(:,name) = S(:,2);
        S_allZ(:,name) = S(:,3);
    end


    %% Plot backbone positions of all samples (averaged)
    
    % First plot this for the figure that will be saved for each individual cluster 
    varMetricS = plotWireframeMeanVar(S_allX,S_allY,S_allZ);
    %spatialVarUpperLim = max(varMetricS);
    a = colorbar; a.Label.String = 'Spatial variance (Angstroms)'; a.Label.FontSize = 15;
    clim([0 spatialVarUpperLim])
    colormap(pink)
    %set(gca,'ColorScale','log')
    set(gcf,'Color',[0.7 0.7 0.7])
    saveas(gca,[folderName,'/outputs/models/cluster',num2str(clust),'_wireframe.fig'])
    close %Close the figure once it's saved
    
    % Next plot this, sans individual wireframes, as a subplot in the output gui 
    nexttile([2 1])
    hold all
    title(['Cluster ',num2str(clust),'; # of structures = ',num2str(nStructures)],'Color',colorsSpaced(clust+1,:))
    varMetricS = plotWireframeMeanVar(S_allX,S_allY,S_allZ);
    clim([0 spatialVarUpperLim])

    clear S_allX S_allY S_allZ %These will change in size depending on number of structures in each cluster
end

    a = colorbar; a.Label.String = 'Spatial variance (Angstroms)'; a.Label.FontSize = 15;
    %clim([0 spatialVarUpperLim])
    colormap(pink)
    %set(gca,'ColorScale','log')

saveas(gcf,[folderName,'/outputs/ClassAveragedEnsemble.fig'])


%% Isolate the mean backbone plot and include sequence information 
%%  Not yet implemented. It's surprisingly annoying to plot letters in specific places...
%
% if ~isempty(seq)
%     seq = num2cell(seq);
%     % Make each base a different color
%     for i = 1:numel(seq)
%         if seq{i} == 'A'
%             %seqColor{i} = [0.9 0.6 0.1];
%             seqColor{i} = [1 0.8 0.05];
%         elseif seq{i} == 'C'
%             seqColor{i} = [0.4 0.5 0.1];
%         elseif seq{i} == 'G'
%             seqColor{i} = [0.9 0.1 0.15];
%         elseif seq{i} == 'U' || seq{i} == 'T'
%             seqColor{i} = [0.1 0.6 0.9];
%         end
%     end
%     lineColor = [0.6 0.6 0.6];
%     font = 'Poor Richard'; % poor, poor Richard :'(
%     %font = 'Courier';
% 
%     % Split off first and last base in sequence
%     seqStart = seq{1};
%     seqMid = {seq{2:numel(seq)-1}};
%     seqEnd = seq{end};
%     seqColorStart = seqColor{1};
%     seqColorMid = {seqColor{2:numel(seq)-1}};
%     seqColorEnd = seqColor{end};
% 
%     figure; hold all
%     set(gcf,'color','w')
%     set(gcf,'Position',[850 50 600 875])
%     %scatter3(meanS(:,1),meanS(:,2),meanS(:,3),100,lineColor,'filled','MarkerEdgeColor','K','LineWidth',1)
%     plot3(meanS(:,1),meanS(:,2),meanS(:,3),'Color',lineColor,'LineWidth',1.5,'LineStyle','--')
% 
%     % Place 5' base a bit before the chain
%     vStart = [meanS(2,3)-meanS(1,3),meanS(2,2)-meanS(1,2),meanS(2,1)-meanS(1,1)];
%     text(meanS(1,1)-vStart(1), meanS(1,2)-vStart(2), meanS(1,3)-vStart(3),...
%         seqStart,'FontSize',35,'Color',seqColorStart,'FontName',font)
% 
%     for i = 1:numel(seqMid)
%         text(mean([meanS(i,1),meanS(i+1,1)]), mean([meanS(i,2),meanS(i+1,2)]), mean([meanS(i,3),meanS(i+1,3)]),...
%             seqMid{i},'FontSize',35,'Color',seqColorMid{i},'FontName',font)
%     end
% 
%     % Place 3' base a bit after the chain
%     vEnd = [meanS(end,3)-meanS(end-1,3),meanS(end,2)-meanS(end-1,2),meanS(end,1)-meanS(end-1,1)];
%     text(meanS(end,1)+vEnd(1), meanS(end,2)+vEnd(2), meanS(end,3)+vEnd(3),...
%         seqEnd,'FontSize',35,'Color',seqColorEnd,'FontName',font)
% 
% end


disp('The class averaging result is saved in .../outputs.')
disp('You may now look at each individual conformer with the MATLAB figure viewer and scale/rotate to your liking.')

end

