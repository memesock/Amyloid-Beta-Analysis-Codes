function alignClusters_allToOne(folderName,files_grouped,nClusters,clusterPDBnumbers)
%% Align PDBs within a spectral cluster, and record all pairwise RMSD values
%   Adapted from alignAll.pml, originally written in Python for use in PyMol
%      MATLAB 2023 has 'pdbsuperpose' function, which should be similar to
%      the Align functionality in PyMol.
%
%   Designed for ensembles of disordered macromolecules
%
%   GW - 2024
%
%   Doing this the way it's done in alignCluster.pml doesn't seem to work,
%   here I'll try aligning all pdbs simply to the first pdb.



%output = fopen([folderName,'/rmsd.txt'], 'wt');



for clust = 0:(nClusters-1)
    files_inThisCluster = files_grouped{clust+1};
    clusterPDBnumbers_thisCluster = clusterPDBnumbers{clust+1};

    % rmsds = zeros([numel(files_inThisCluster),numel(files_inThisCluster)]);
    % 
    % for i = 1:numel(files_inThisCluster)
    %     %print(["alignment template: structure",num2str(i)])
    % 
    %     for j = 1:numel(files_inThisCluster)
    % 
    %         if i == j % If same structure, RMSD=0
    %             %fprintf(output, '%s\n', [num2str(i),'0',num2str(j),' 0']);
    %             rmsds(i,j) = 0;
    % 
    %         else
    %             pdb1 = pdbread([folderName,'/',files_inThisCluster{i}]);
    %             pdb2 = pdbread([folderName,'/',files_inThisCluster{j}]);
    % 
    %             if isfield(pdb1,'Sequence') == 1 % PDB likely a traditional protein pdb
    %                 [DIST,RMSD] = pdbsuperpose(pdb1,pdb2);
    %             else % PDB likely not a protein pdb (eg nucleic acid, artificial construct)
    %                 [DIST,RMSD] = pdbsuperpose_general(pdb1,pdb2,'SEQALIGN', 'False');
    % 
    %             end
    % 
    %             rmsds(i,j) = RMSD;
    %             %fprintf(output, '%s\n', [num2str(i),'0',num2str(j),' ',num2str(RMSD)]);
    % 
    %         end
    %     end
    % end
    % 
    % rmsdMatrix = rmsds;
    % 
    % 
    % %% Determine minimum RMSD values
    % whereToAlign = zeros([numel(files_inThisCluster),2]);
    % whereToAlign(:,1) = 1:1:numel(files_inThisCluster);
    % 
    % for k = 1:numel(files_inThisCluster)
    % 	currentRow = rmsdMatrix(k,:);
    % 	currentRow(currentRow==0) = NaN;
    % 	minValue = min(currentRow);
    % 	minIndices = find(currentRow==minValue);
    % 	whereToAlign(k,2) = minIndices(1);
    % end

    for i = 1:numel(files_inThisCluster)
    %% Align structures that have minimum resulting RMSDs
        %structure1 = whereToAlign(i,1);
        structure2 = i;

        pdb1 = pdbread([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/',files_inThisCluster{1}]);
        pdb2 = pdbread([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/',files_inThisCluster{structure2}]);

        if isfield(pdb1,'Sequence') == 1 % PDB likely a traditional protein pdb
            [DIST,RMSD,TRANSF,PDB2TX] = pdbsuperpose(pdb1,pdb2);
        else % PDB likely not a protein pdb (eg nucleic acid, coarse grained)
            [DIST,RMSD,TRANSF,PDB2TX] = pdbsuperpose_general(pdb1,pdb2,'SEQALIGN', 'False');
        end

        %data = cmd.align( basename+str(PDBnumbers[structure1]), basename+str(PDBnumbers[structure2]) );
    
        % Save aligned PDBs
        if ~exist([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/Aligned'],'dir')
            mkdir([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/Aligned'])
        end

        %pdbwrite([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/Aligned/',files_inThisCluster{structure2}], PDB2TX)
        
        % If the PDB2TX structure has other shit in it other than 'Model', it can fuck up pdbwrite for some
        % reason. So force it just to write a structure with 'Model' only. 
        saveme.Model = PDB2TX.Model;
        pdbwrite([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/Aligned/',files_inThisCluster{structure2}], saveme)
        %clear DIST RMSD TRANSF PDB2TX

    end

    % for i = 1:numel(files_inThisCluster)
    %    pdbwrite([folderName,'/PDBs_SpectralClustered/Cluster',num2str(clust),'/',files_inThisCluster{i}], PDB2TX)
    % end

end


end

