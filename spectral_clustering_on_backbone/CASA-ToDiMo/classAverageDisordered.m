clear
close all

%% Assign structures within an ensemble to unique classes based on nearest-neighbor RMSD minimization
%   This workflow requires a set of .pdb files in a subfolder of the
%   working directory
%
%   GW - started March 2023
%      - January 2024: moved all Pymol scripting into MATLAB so you don't
%      need to alternate between the programs
%
%      - February 2024: added automatic determination of rmsdThreshold and
%      # of K-means clusters
%
%      - June 2024: Automated the program, and consolidated all plots into
%      a couple of windows. Results are also automatically saved. 
%
%   Ordered list of plots:
%       0) The spectral clustering result (K-means clustered graph) 
%       1) Pairwise RMSD heatmap
%       2) Binary pairwise RMSD matrix (adjacency matrix)
%       3) Entropy metric
%       4) Graph colored by closeness centrality
%       5) Graph colored by degree centrality
%   If rmsdThreshold = 0 is specified and the program determines an
%   RMSDthreshold automatically, two additional plots will be output
%   assessing quality across a spectrum of RMSDthreshold values
%
%


%% User-input variables
% Variables that commonly need to be changed
folderName = 'pdb_10_stride_com_2.5'
% ExamplePool_MixedSequenceSingleStrandedRNA';

rmsdThreshold = 0; % Threshold RMSD value to be considered in the same class, in Angstroms; set to 0 for program to determine for you via Okayness maximization
nKmeans = 0; % Number of clusters to use for final K-means on the node/edge graph; set to 0 and the program will do a search between N_clusters=2-20 and choose N_clusters that maximizes the average silhouette

species = 'protein'; % Currently supported: either 'RNA', 'PAR', 'protein', or 'CG' (already coarse grained)

spatialVarUpperLim = 100; % how high to set the colorbar when visualizing the classed averaged conformers; 100A by default

seq = []; % NOTE: sequence functionality has not been implemented yet. But when it is: 
% If you want sequence info to be included when plotting conformers, enter the sequence; eg:
%      "sequence = ['SRKMCTVL']" for an amino acid sequence
%      "sequence = ['AUGCUCG']" for an RNA sequence
%      Leave "sequence = []" if you don't want to plot with sequence info present


PARtrueOrder = [1;2;13;16;17;18;19;20;21;22;3;4;5;6;7;8;9;10;11;12;14;15];
%PARtrueOrder = [1;2;8;9;10;11;12;13;14;15;3;4;5;6;7]; % If species='PAR', need to manually look in the pdb and enter correct subunit order (For some reason the PAR subunit indices get reordered);
% if species is not 'PAR', this variable does not matter.

rmsdThresholdSearchRange = (1.00:0.2:24.0); % range of RMSD values that will be tried if rmsdThreshold=0 (in Anstroms); this does not usually need to be changed




%%
disp('----------------------------------------------------------------------------------------')
disp('CASA-ToDiMo: Class Averaging via Spectral Analysis of Totally Disordered Macromolecules')
disp('Cornell University')
disp('----------------------------------------------------------------------------------------')

%% Define subplot layout
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(4,4)

%% Perform initial alignment and pairwise RMSD calculations

disp(' ')
if exist([folderName,'/','rmsdMatrix.mat'],'file') ~= 0 % Don't do again if RMSDs already computed
    disp('Loading previously computed pairwise RMSDs...')
    load([folderName,'/','rmsdMatrix.mat'])
    Nstructures = sqrt(numel(rmsdMatrix));
    %rmsdMatrix = reshape(rmsdMatrix,[Nstructures, Nstructures]);
else
    disp('Computing RMSDs between each pair of structures after alignment... This may take a while.')

    packages = ver().Name; % Check if Parallel toolbox is installed; if so do it the faster way
    if max(strcmpi(packages,'Parallel Computing Toolbox')) == 1
        rmsdMatrix = alignAll_parallel(folderName);
    else
        rmsdMatrix = alignAll(folderName);
    end

    save([folderName,'/rmsdMatrix'],'rmsdMatrix') % save for future reruns
    Nstructures = sqrt(numel(rmsdMatrix));
end
disp('Done.')


%% Automatically determine RMSD threshold if desired, based on okayness metric
if rmsdThreshold == 0
    disp(['[', datestr(now,'HH:MM:SS'), '] Entering determineRMSDthresh'])
    rmsdThreshold = determineRMSDthresh(Nstructures, rmsdThresholdSearchRange, rmsdMatrix);
    disp(['[', datestr(now,'HH:MM:SS'), '] Finished determineRMSDthresh'])
end

%% Define binary RMSD pairings (connectivity matrix)
RMSD_binary = zeros(Nstructures, Nstructures);
for i = 1:Nstructures
    for j = 1:Nstructures
        if rmsdMatrix(i,j) <= rmsdThreshold
            RMSD_binary(i,j) = 1;
        else
            RMSD_binary(i,j) = 0;
        end
    end
end


%% Perform spectral clustering on connectivity matrix; plot graph

C = SpecClust(RMSD_binary, nKmeans);
clusters = [1:Nstructures; C']';
clusters = sortrows(clusters,2);

% Color graph by K-means clusters identified
%figure('Name','Clustered Graph'); hold all
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
nClusters = numel(unique(clusters(:,2)));
colors = colormap(jet);
colorSpacing = floor(numel(colors(:,1)) / (nClusters+1) * 0.9);
for j = 1:nClusters+1
    colorsSpaced(j,:) = colors(j*colorSpacing,:);

    clusterIndices{j} = clusters(clusters(:,2)==j-1);
    highlight(h,clusterIndices{j},'NodeColor',colorsSpaced(j,:))
end


%% Optional plots (specified by plotsToInclude):
%% Plot pairwise RMSD results

% Plot heatmap
%figure('Renderer', 'painters', 'Position', [10 10 580 540])
nexttile([2 1])
h = heatmap(rmsdMatrix);
h.GridVisible = 'off';
colormap('jet')
set(gcf,'color','w')
set(gca,'FontSize',12)
Labels = 1:Nstructures; CustomLabels = string(Labels);
CustomLabels(mod(Labels,Nstructures) ~= 0) = " "; % Replace all but the last elements with spaces, for clarity
h.XDisplayLabels = CustomLabels; h.YDisplayLabels = CustomLabels;
title('RMSD between all pairs of structures')
%xlabel('structure'); ylabel('structure')
%annotation('textarrow',[1,1],[0.5,0.5],'string','RMSD','Interpreter','latex','Fontsize',20, ...
% 'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');



% Plot connectivity matrix
%figure('Renderer', 'painters', 'Position', [10 10 580 540])
nexttile([2 1])
h = heatmap(RMSD_binary);
h.GridVisible = 'off';
colormap('jet')
%clim([0,1])
set(gcf,'color','w')
set(gca,'FontSize',12)
h.XDisplayLabels = CustomLabels; h.YDisplayLabels = CustomLabels;
h.ColorbarVisible = 'off';
title('Binarized RMSD matrix (adjacency matrix)')
%xlabel('structure'); ylabel('structure')



%% Plot disorder metrics
numRed = numel(find(RMSD_binary==1)) - Nstructures; % subtract diagonals (same structures)
fracRed = numRed/(Nstructures*Nstructures);
numBlue = numel(find(RMSD_binary==0))+1; % add one starting state
fracBlue = numBlue/(Nstructures*Nstructures);
for i = 1:Nstructures
    histRed(i) = numel(find(RMSD_binary(:,i)==1));
    histBlue(i) = numel(find(RMSD_binary(:,i)==0));
end
histRed = histRed - 1; % subtract diagonals (same structures)
histBlue = histBlue + 1; % add one starting state
fracHistRed = histRed ./ Nstructures;
fracHistBlue = histBlue ./ Nstructures; 

kB = 1.380649E-23; % Boltzmann constant, in J/K
W = floor(numBlue/2); % Assume two structures occupy different microstates if they are disconnected
Sdivk = log(W);

%figure('Renderer', 'painters', 'Position', [10 10 580 540])
nexttile([2 1])
box on
set(gcf,'color','w')
set(gca,'FontSize',20,'LineWidth',2)
%histogram(fracHistRed,10,'FaceColor',[0.6350 0.0780 0.1840],'LineWidth',2,'FaceAlpha',1.0)
histogram(fracHistBlue,10,'FaceColor','b','LineWidth',2,'FaceAlpha',1.0)
xlabel('Fraction of disconnected structures')
ylabel('# of structures in ensemble')
xlim([0.5 1])
title(['Total fraction of disconnected structures = ',num2str(fracBlue),'; S/$k_{B}$ = ',num2str(Sdivk)],'FontSize',12,...
    'Interpreter','latex')


%% Plot graph colored by closeness centrality
%figure('Name','Centrality Graph'); hold all
nexttile([2 1])
G = graph(RMSD_binary,'omitselfloops');
cent = centrality(G,'closeness');

h = plot(G,'-db','LineWidth',1,'MarkerSize',5);
set(gcf,'color','w')
set(gca,'FontSize',12)
set(gca,'LineWidth',2)
box on
grid off
h.EdgeColor = [0.5 0.5 0.5];
h.Marker = 'O';

h.NodeCData = cent;
colormap jet
a = colorbar;
ylabel(a,'Closeness Centrality','FontSize',16,'Rotation',270);
clim([0,6E-3])

title(['Mean centrality = ',num2str(mean(cent))])


%% Plot graph colored by degree centrality
%figure('Name','Centrality Graph'); hold all
nexttile([2 1])
G = graph(RMSD_binary,'omitselfloops');
cent = centrality(G,'degree');

h = plot(G,'-db','LineWidth',1,'MarkerSize',5);
set(gcf,'color','w')
set(gca,'FontSize',12)
set(gca,'LineWidth',2)
box on
grid off
h.EdgeColor = [0.5 0.5 0.5];
h.Marker = 'O';

h.NodeCData = cent;
colormap jet
a = colorbar;
ylabel(a,'Degree Centrality','FontSize',16,'Rotation',270);
clim([0,35])

title(['Mean centrality = ',num2str(mean(cent))])


%% Assess clustering result and save if satisfied

prompt = "\nSpectral clustering results vary depending on the RMSD threshold and # of K-means clusters declared." + ...
    "\nAre you satisfied with the results? Enter 1 for yes, 0 for no (redo spectral clustering)." + ...
    "\nYou do not need to recompute the pairwise RMSD matrix, so subsequent runs will take much less computing time." + ...
    "\nNote that if many stray points are present in the graph that are completely disconnected from other points, the next step might not work." + ...
    "\n";
%satisfied = input(prompt);
satisfied = 1
if satisfied == 1
    disp('Nice! Saving outputs...')
    if not(isfolder([folderName,'/outputs']))
        mkdir([folderName,'/outputs'])
    end
    if not(isfolder([folderName,'/outputs/models']))
        mkdir([folderName,'/outputs/models'])
    end
    save([folderName,'/outputs/KmeansClusters.mat'],'clusters')
    saveas(gcf,[folderName,'/outputs/SpectralAnalysis.png'])

    disp('Writing to log...')
    writelines(['3D class averaging on structural ensemble in folder: ',folderName],[folderName,'/outputs/log.txt'])
    writelines(['Species: ',species],[folderName,'/outputs/log.txt'],WriteMode='append')
    writelines('Pairwise RMSD values saved in rmsd.txt',[folderName,'/outputs/log.txt'],WriteMode='append')
    writelines(['Number of conformers: ',num2str(Nstructures)],[folderName,'/outputs/log.txt'],WriteMode='append')
    writelines(['RMSD threshold for graph calculation: ',num2str(rmsdThreshold)],[folderName,'/outputs/log.txt'],WriteMode='append')
    writelines(['Number of clusters: ',num2str(nClusters)],[folderName,'/outputs/log.txt'],WriteMode='append')

    % Plot class averaged conformers (a new figure window will open)
    disp('The program will now compute the class averaged conformers and display them when done...')
    disp('----------------------------------------------------------------------------------------')

    SphagettiPlot(folderName, species, seq, RMSD_binary, G, clusters, spatialVarUpperLim, PARtrueOrder)

elseif satisfied == 0
    error('Feel free to change the rmsdThreshold and/or nKmeans and run "classAverageDisordered.m" again.')

end


