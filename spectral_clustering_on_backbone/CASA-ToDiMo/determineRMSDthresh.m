function rmsd_threshold = determineRMSDthresh(Nstructures, rmsdThresholdSearchRange,RMSD)

%% Assign structures within an ensemble to unique classes based on nearest-neighbor RMSD minimization 
%   This workflow requires you to first generate RMSD values of all possible pairwise alignments
%       using PYMOL (can use 'alignAll.pml' script)
%
%   Designed for ensembles of disordered nucleic acid conformers
%   
%  Determine the optimal RMSDthreshold value toward generating an adjacency
%  matrix from a pairwise RMSD matrix by iteratively optimizing an
%  'Okayness' metric 


%% Load RMSDs and set up pairwise RMSD matrix 
%RMSD_load = readmatrix([folderName,'/rmsd.txt']);
%RMSD = RMSD_load(:,2);
%RMSD = reshape(RMSD,[Nstructures, Nstructures]);


for k = 1:numel(rmsdThresholdSearchRange)

    rmsdThreshold = rmsdThresholdSearchRange(k);
    for i = 1:Nstructures
        for j = 1:Nstructures
            if RMSD(i,j) <= rmsdThreshold
                RMSD_binary(i,j) = 1;
            else
                RMSD_binary(i,j) = 0;
            end
        end
    end


    %% Define graph and compute Okayness
  
    G = graph(RMSD_binary,'omitselfloops','upper');
    closenessCent = centrality(G,'closeness');

    nEdges(k) = numel(G.Edges)/2;
    edgeConnectivity(k) = computeEdgeConn(G);
    %okayness(k) = 1/(sqrt(nEdges(k))).*mean(closenessCent);
    %okayness(k) = 1/(sqrt(nEdges(k))).*mean(closenessCent)./(1+0.1*edgeConnectivity(k));
    okayness(k) = 1/(sqrt(nEdges(k))).*mean(closenessCent)./(1+log(1+edgeConnectivity(k)));

end

%figure; hold all
nexttile([2 1])
box on; grid on
set(gcf,'color','w')
set(gca,'LineWidth',2,'FontSize',15)
plot(rmsdThresholdSearchRange,nEdges,'ko-','LineWidth',2)
yyaxis right 
%plot(rmsdThresholdSearchRange,edgeConnectivity,'ro-','LineWidth',2)
xlim([rmsdThresholdSearchRange(1),rmsdThresholdSearchRange(end)])
xlabel('RMSD threshold (\AA)','Interpreter','latex','FontSize',22)
legend('# of edges','edge connectivity','Location','NorthWest')

%figure; hold all
nexttile([2 1])
box on; grid on
set(gcf,'color','w')
set(gca,'LineWidth',2,'FontSize',15)
plot(rmsdThresholdSearchRange,okayness,'o-','LineWidth',2)
xlim([rmsdThresholdSearchRange(1),rmsdThresholdSearchRange(end)])
xlabel('RMSD threshold (\AA)','Interpreter','latex','FontSize',22)
ylabel('Okayness','FontSize',22)


%% Determine RMSDthreshold corresponding to maximal okayness
[okayness_max, whereMax] = max(okayness);
rmsd_threshold = rmsdThresholdSearchRange(whereMax)

end


%% Auxiliary functions
function k = computeEdgeConn(g) % From Christine Tobler (Mathworks)
    % Make sure the graph is unweighted:
    if contains("Weight", g.Edges.Properties.VariableNames)
        g.Edges.Weight = [];
    end
    % Compute the maximum flow from node 1 to every other node
    mf = zeros(1, numnodes(g)-1);
    for ii=2:numnodes(g)
        mf(ii-1) = maxflow(g, 1, ii);
    end
    % The edge connectivity is the minimum of these maximum flows
    k = min(mf);
end
