function [files_grouped, nClusters, clusterPDBnumbers] = GroupClusters(folderName, clusters)

%% Copy all PDBs into subfolders in the working directory based on cluster number 
%   Only execute this part when you are satisfied with the clustering situation 
%
%   Requires you to input the 'clusters' variable [Nstructuresx2 matrix]
%   from 'RNA_classAverage.m' script and load it into the workspace before
%   running. 
%
%   GW - 2023 March  


% basename = 'frame'; % preceeding the number 
% postname = '.pdb'; % following the number 

files = dir( fullfile(folderName,'*.pdb') );
files = {files.name}';
%files = natsort(files); % Sort filenames in numerical order; this line requires an external script so I just disabled it, it's mostly just organizational


%% Do the file copying 
nClusters = numel(unique(clusters(:,2)));
clusterPDBsAll = clusters(:,1);

if ~exist([folderName,'/PDBs_SpectralClustered'],'dir')
    mkdir([folderName,'/PDBs_SpectralClustered'])
end

for i = 0:(nClusters-1)
    if ~exist([folderName,'/PDBs_SpectralClustered/Cluster',num2str(i)],'dir')
        mkdir([folderName,'/PDBs_SpectralClustered/Cluster',num2str(i)])
    end
    clusterPDBindices = find(clusters(:,2)==i);
    clusterPDBnumbers{i+1} = clusterPDBsAll(clusterPDBindices);
    
    clusterPDBnumbers_thisCluster = clusterPDBnumbers{i+1};
    files_grouped{i+1} = {files{clusterPDBnumbers_thisCluster}}';

    for j = 1:numel(clusterPDBnumbers_thisCluster)
        copyfile([folderName,'/',files{clusterPDBnumbers_thisCluster(j)}],...
            [folderName,'/PDBs_SpectralClustered/Cluster',num2str(i)])
    end

end



end








