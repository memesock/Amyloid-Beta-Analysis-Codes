function [S] = getWireframeIndividual(pdb,species,PARtrueOrder)
%Coarse grain macromolecules by sampling positions of key atoms.

if ~strcmpi(species,'CG') % get names of atoms if the model is not CG
    for i = 1:pdb(1,end).AtomSerNo
        atomNames{i} = pdb(1,i).AtomName;
    end
end

%%

if strcmpi(species,'RNA') %For RNA, sample phosphorus atoms along backbone
    % Get coordinates of phosphorus atoms
    wherePhos = find(strcmpi(atomNames,'P'));
    for i = 1:numel(wherePhos)
        Ps_x(i,1) = pdb(1,wherePhos(i)).X ;
        Ps_y(i,1) = pdb(1,wherePhos(i)).Y ;
        Ps_z(i,1) = pdb(1,wherePhos(i)).Z ;
    end

    S(:,1) = Ps_x;
    S(:,2) = Ps_y;
    S(:,3) = Ps_z;

elseif strcmpi(species,'PAR') % Use poly(ADP-ribose) sampling method
    % Sample average of each pair of P atoms
    wherePhos = [find(strcmpi(atomNames,'PA')), find(strcmpi(atomNames,'PB'))]; wherePhos = sort(wherePhos);
    for i = 1:2:numel(wherePhos)
        Ps_x(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).X, pdb(1,wherePhos(i+1)).X] ) ;
        Ps_y(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).Y, pdb(1,wherePhos(i+1)).Y] ) ;
        Ps_z(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).Z, pdb(1,wherePhos(i+1)).Z] ) ;
    end

    %Reorder
    Ps_x_sorted = Ps_x(PARtrueOrder); Ps_y_sorted = Ps_y(PARtrueOrder); Ps_z_sorted = Ps_z(PARtrueOrder);

    S(:,1) = Ps_x_sorted;
    S(:,2) = Ps_y_sorted;
    S(:,3) = Ps_z_sorted;

elseif strcmpi(species,'protein') % For proteins/peptides, sample alpha carbons along backbone
    % Get coordinates of alpha carbons
    whereAlph = find(strcmpi(atomNames,'BB'));
    for i = 1:numel(whereAlph)
        Ps_x(i,1) = pdb(1,whereAlph(i)).X ;
        Ps_y(i,1) = pdb(1,whereAlph(i)).Y ;
        Ps_z(i,1) = pdb(1,whereAlph(i)).Z ;
    end

    S(:,1) = Ps_x;
    S(:,2) = Ps_y;
    S(:,3) = Ps_z;

elseif strcmpi(species,'CG') % Models already coarse grained, simply sample every atom
    % for i = 1:numel(pdb)
    %     atoms_x(i,1) = pdb.X;
    %     atoms_y(i,1) = pdb.Y;
    %     atoms_z(i,1) = pdb.Z;
    % end
    % 
    % S(:,1) = atoms_x;
    % S(:,2) = atoms_y;
    % S(:,3) = atoms_z;

    S(:,1) = [pdb.X]';
    S(:,2) = [pdb.Y]';
    S(:,3) = [pdb.Z]';
end



end
