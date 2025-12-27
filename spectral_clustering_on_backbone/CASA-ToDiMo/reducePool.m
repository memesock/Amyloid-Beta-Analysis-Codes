clear; close all
%% Copy every a'th pdb file for even pool reduction
%
%   GW - September 2024
%

poolFolder = 'rU30_SAXSMD_HBCUFIX_2000PDBs_50ns';
outputFolder = 'rU30_SAXSMD_HBCUFIX_100PDBs_50ns';
nFinal = 100; % How many pdb files you want in the final pool


%% 
pool = dir([poolFolder,'/*.pdb']);
nPool = numel(pool);
a = nPool/nFinal;

% Sort pool names in ascending numerical order (Cedric from Matlab help)
[~, reindex] = sort( str2double( regexp( {pool.name}, '\d+', 'match', 'once' )));
pool = pool(reindex) ;


poolFiles = {pool.name};
reducedPoolFiles = poolFiles(1:a:nPool);

if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end

for i = 1:nFinal
    copyfile([poolFolder,'/',reducedPoolFiles{i}],[outputFolder,'/',reducedPoolFiles{i}])
end



