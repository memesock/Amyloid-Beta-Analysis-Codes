clear
close all

t_TOTAL = tic;
fprintf('\n[%s] Starting CASA-ToDiMo pipeline\n', datestr(now,'HH:MM:SS'));

%% Assign structures within an ensemble to unique classes based on nearest-neighbor RMSD minimization
% (comments unchanged)

%% User-input variables
folderName = 'pdb_10_stride_com_2.5'

rmsdThreshold = 0;
nKmeans = 0;

species = 'protein';
spatialVarUpperLim = 100;
seq = [];

PARtrueOrder = [1;2;13;16;17;18;19;20;21;22;3;4;5;6;7;8;9;10;11;12;14;15];

rmsdThresholdSearchRange = (5.0:1.0:25.0);

disp('----------------------------------------------------------------------------------------')
disp('CASA-ToDiMo: Class Averaging via Spectral Analysis of Totally Disordered Macromolecules')
disp('Cornell University')
disp('----------------------------------------------------------------------------------------')

%% Define subplot layout
fprintf('[%s] Initializing figure layout\n', datestr(now,'HH:MM:SS'));
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(4,4)

%% Perform initial alignment and pairwise RMSD calculations
disp(' ')
t_RMSD = tic;

if exist([folderName,'/','rmsdMatrix.mat'],'file') ~= 0
    disp('Loading previously computed pairwise RMSDs...')
    load([folderName,'/','rmsdMatrix.mat'])
    Nstructures = sqrt(numel(rmsdMatrix));
else
    disp('Computing RMSDs between each pair of structures after alignment... This may take a while.')

    packages = ver().Name;
    if max(strcmpi(packages,'Parallel Computing Toolbox')) == 1
        fprintf('[%s] Using parallel RMSD alignment\n', datestr(now,'HH:MM:SS'));
        rmsdMatrix = alignAll_parallel(folderName);
    else
        fprintf('[%s] Using serial RMSD alignment\n', datestr(now,'HH:MM:SS'));
        rmsdMatrix = alignAll(folderName);
    end

    save([folderName,'/rmsdMatrix'],'rmsdMatrix')
    Nstructures = sqrt(numel(rmsdMatrix));
end

fprintf('[%s] RMSD matrix ready | Nstructures = %d | Time = %.2f min\n', ...
    datestr(now,'HH:MM:SS'), Nstructures, toc(t_RMSD)/60);
disp('Done.')

%% Automatically determine RMSD threshold
if rmsdThreshold == 0
    fprintf('[%s] Entering determineRMSDthresh\n', datestr(now,'HH:MM:SS'));
    t_thresh = tic;
    rmsdThreshold = determineRMSDthresh(Nstructures, rmsdThresholdSearchRange, rmsdMatrix);
    fprintf('[%s] Finished determineRMSDthresh | RMSD threshold = %.3f | Time = %.2f min\n', ...
        datestr(now,'HH:MM:SS'), rmsdThreshold, toc(t_thresh)/60);
end

%% Define binary RMSD pairings
fprintf('[%s] Building binary RMSD adjacency matrix\n', datestr(now,'HH:MM:SS'));
t_bin = tic;

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

fprintf('[%s] Adjacency matrix built | Time = %.2f s\n', ...
    datestr(now,'HH:MM:SS'), toc(t_bin));

%% Perform spectral clustering
fprintf('[%s] Entering spectral clustering\n', datestr(now,'HH:MM:SS'));
t_spec = tic;

C = SpecClust(RMSD_binary, nKmeans);
clusters = [1:Nstructures; C']';
clusters = sortrows(clusters,2);

fprintf('[%s] Spectral clustering finished | Time = %.2f min\n', ...
    datestr(now,'HH:MM:SS'), toc(t_spec)/60);

%% Plot clustered graph
fprintf('[%s] Plotting clustered graph\n', datestr(now,'HH:MM:SS'));

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

nClusters = numel(unique(clusters(:,2)));
colors = colormap(jet);
colorSpacing = floor(numel(colors(:,1)) / (nClusters+1) * 0.9);

for j = 1:nClusters+1
    colorsSpaced(j,:) = colors(j*colorSpacing,:);
    clusterIndices{j} = clusters(clusters(:,2)==j-1);
    highlight(h,clusterIndices{j},'NodeColor',colorsSpaced(j,:))
end

%% Plot RMSD heatmap
fprintf('[%s] Plotting RMSD heatmap\n', datestr(now,'HH:MM:SS'));
nexttile([2 1])
h = heatmap(rmsdMatrix);
h.GridVisible = 'off';
colormap('jet')
set(gcf,'color','w')
set(gca,'FontSize',12)

Labels = 1:Nstructures;
CustomLabels = string(Labels);
CustomLabels(mod(Labels,Nstructures) ~= 0) = " ";
h.XDisplayLabels = CustomLabels;
h.YDisplayLabels = CustomLabels;
title('RMSD between all pairs of structures')

%% Plot connectivity matrix
fprintf('[%s] Plotting adjacency matrix\n', datestr(now,'HH:MM:SS'));
nexttile([2 1])
h = heatmap(RMSD_binary);
h.GridVisible = 'off';
colormap('jet')
set(gcf,'color','w')
set(gca,'FontSize',12)
h.XDisplayLabels = CustomLabels;
h.YDisplayLabels = CustomLabels;
h.ColorbarVisible = 'off';
title('Binarized RMSD matrix (adjacency matrix)')

%% Disorder metrics
fprintf('[%s] Computing disorder metrics\n', datestr(now,'HH:MM:SS'));
t_dis = tic;

numRed = numel(find(RMSD_binary==1)) - Nstructures;
fracRed = numRed/(Nstructures*Nstructures);
numBlue = numel(find(RMSD_binary==0))+1;
fracBlue = numBlue/(Nstructures*Nstructures);

for i = 1:Nstructures
    histRed(i) = numel(find(RMSD_binary(:,i)==1));
    histBlue(i) = numel(find(RMSD_binary(:,i)==0));
end

histRed = histRed - 1;
histBlue = histBlue + 1;

fracHistRed = histRed ./ Nstructures;
fracHistBlue = histBlue ./ Nstructures;

W = floor(numBlue/2);
Sdivk = log(W);

fprintf('[%s] Disorder metrics done | Time = %.2f s\n', ...
    datestr(now,'HH:MM:SS'), toc(t_dis));

nexttile([2 1])
box on
set(gcf,'color','w')
set(gca,'FontSize',20,'LineWidth',2)
histogram(fracHistBlue,10,'FaceColor','b','LineWidth',2,'FaceAlpha',1.0)
xlabel('Fraction of disconnected structures')
ylabel('# of structures in ensemble')
xlim([0.5 1])
title(['Total fraction of disconnected structures = ',num2str(fracBlue),'; S/$k_{B}$ = ',num2str(Sdivk)], ...
    'FontSize',12,'Interpreter','latex')

%% Closeness centrality graph
fprintf('[%s] Computing closeness centrality\n', datestr(now,'HH:MM:SS'));
t_cent = tic;

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

fprintf('[%s] Closeness centrality done | Time = %.2f s\n', ...
    datestr(now,'HH:MM:SS'), toc(t_cent));

%% Degree centrality graph
fprintf('[%s] Computing degree centrality\n', datestr(now,'HH:MM:SS'));
t_deg = tic;

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

fprintf('[%s] Degree centrality done | Time = %.2f s\n', ...
    datestr(now,'HH:MM:SS'), toc(t_deg));

%% Save results
disp('Nice! Saving outputs...')
t_save = tic;

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
writelines(['Number of conformers: ',num2str(Nstructures)],[folderName,'/outputs/log.txt'],WriteMode='append')
writelines(['RMSD threshold: ',num2str(rmsdThreshold)],[folderName,'/outputs/log.txt'],WriteMode='append')
writelines(['Number of clusters: ',num2str(nClusters)],[folderName,'/outputs/log.txt'],WriteMode='append')

fprintf('[%s] Saving done | Time = %.2f s\n', ...
    datestr(now,'HH:MM:SS'), toc(t_save));

disp('The program will now compute the class averaged conformers and display them when done...')
disp('----------------------------------------------------------------------------------------')

t_spag = tic;
SphagettiPlot(folderName, species, seq, RMSD_binary, G, clusters, spatialVarUpperLim, PARtrueOrder)
fprintf('[%s] SphagettiPlot finished | Time = %.2f min\n', ...
    datestr(now,'HH:MM:SS'), toc(t_spag)/60);

fprintf('\n[%s] TOTAL WALL TIME = %.2f minutes\n', ...
    datestr(now,'HH:MM:SS'), toc(t_TOTAL)/60);

