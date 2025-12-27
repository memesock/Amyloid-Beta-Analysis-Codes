clear
close all

%% Assign structures within an ensemble to unique classes based on nearest-neighbor RMSD minimization
%   This script just computes the rmsd matrices for a bunch of files in a row, so you can do this overnite and do the 
%       rest of the analysis later
%
%   GW - September 2024 
%


%% User-input variables
% Variables that commonly need to be changed
folderName = {'rC30_SAXSMD_HBCUFIX_100PDBs_60ns','rU30_SAXSMD_HBCUFIX_100PDBs_50ns','rC30_SAXSMD_HBCUFIX_500PDBs_60ns'};

rmsdThreshold = 0; % Threshold RMSD value to be considered in the same class, in Angstroms; set to 0 for program to determine for you via Okayness maximization
%nKmeans = 4; % Number of clusters to use for final K-means on the node/edge graph; set to 0 and the program will do a search between N_clusters=2-20 and choose N_clusters that maximizes the average silhouette

rmsdThresholdSearchRange = (1:0.1:20); % range of RMSD values that will be tried if rmsdThreshold=0 (in Anstroms); this does not usually need to be changed




%%
disp('----------------------------------------------------------------------------------------')
disp('CASA-ToDiMo: Class Averaging via Spectral Analysis of Totally Disordered Macromolecules')
disp('Cornell University')
disp('----------------------------------------------------------------------------------------')



%% Perform initial alignment and pairwise RMSD calculations

disp(' ')

for n = 1:numel(folderName)

    folder = folderName{n}

    % if exist([folder,'/','rmsdMatrix.mat'],'file') ~= 0 % Don't do again if RMSDs already computed
    %     disp('Loading previously computed pairwise RMSDs...')
    %     %RMSD_load = readmatrix([folderName,'/rmsd.txt']);
    %     load([folder,'/','rmsdMatrix.mat'])
    %     %rmsdMatrix = RMSD_load(:,2);
    %     Nstructures = sqrt(numel(rmsdMatrix));
    %     %rmsdMatrix = reshape(rmsdMatrix,[Nstructures, Nstructures]);
    % else
        disp('Computing RMSDs between each pair of structures after alignment... This may take a while.')
        rmsdMatrix = alignAll_parallel(folder);
        save([folder,'/rmsdMatrix'],'rmsdMatrix') % save for future reruns
        Nstructures = sqrt(numel(rmsdMatrix));
    % end
    disp('Done.')


    %% Automatically determine RMSD threshold if desired, based on okayness metric
    if rmsdThreshold == 0
        rmsdThreshold = determineRMSDthresh(Nstructures,rmsdThresholdSearchRange,rmsdMatrix);
    end

end

