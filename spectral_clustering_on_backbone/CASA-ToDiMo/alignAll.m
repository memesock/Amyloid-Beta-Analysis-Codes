function [rmsdMatrix] = alignAll(folderName)
%% Align every pair of PDBs in a directory, and record all pairwise RMSD values
%   Adapted from alignAll.pml, originally written in Python for use in PyMol
%      MATLAB 2023 has 'pdbsuperpose' function, which should be similar to
%      the Align functionality in PyMol.
%
%   Designed for ensembles of disordered macromolecules
%
%   GW - 2024


files = dir( fullfile(folderName,'*.pdb') );
files = {files.name}';
%files = natsort(files); % Sort filenames in numerical order
output = fopen([folderName,'/rmsd.txt'], 'wt');
rmsds = zeros([numel(files),numel(files)]);

% Read pdbs first before executing main loop 
for i = 1:numel(files)
     pdbArray{i} = pdbread([folderName,'/',files{i}]);
end

for i = 1:numel(files)
    %print(["alignment template: structure",num2str(i)])

    for j = 1:numel(files)

        if i == j % If same structure, RMSD=0
            fprintf(output, '%s\n', [num2str(i),'0',num2str(j),' 0']);
            rmsds(i,j) = 0;

        else
            pdb1 = pdbArray{i};
            pdb2 = pdbArray{j}; 

            if isfield(pdb1,'Sequence') == 1 % PDB likely a traditional protein pdb
                [DIST,RMSD,TRANSF,PDB2TX] = pdbsuperpose(pdb1,pdb2);
            else % PDB likely not a protein pdb (eg nucleic acid, coarse grained)
                [DIST,RMSD,TRANSF,PDB2TX] = pdbsuperpose_general(pdb1,pdb2,'SEQALIGN', 'False');

            end

            rmsds(i,j) = RMSD;
            fprintf(output, '%s\n', [num2str(i),'0',num2str(j),' ',num2str(RMSD)]);

        end
    end
end

rmsdMatrix = rmsds;

end

