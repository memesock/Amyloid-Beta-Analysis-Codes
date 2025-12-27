function rmsd_threshold = determineRMSDthresh(Nstructures, rmsdThresholdSearchRange, RMSD)

fprintf('[%s] Entering determineRMSDthresh\n', datestr(now,'HH:MM:SS'));
fprintf('  Nstructures              = %d\n', Nstructures);
fprintf('  # RMSD thresholds to try = %d\n', numel(rmsdThresholdSearchRange));
fprintf('------------------------------------------------------------\n');

t_total = tic;

nT = numel(rmsdThresholdSearchRange);

for k = 1:nT

    rmsdThreshold = rmsdThresholdSearchRange(k);

    fprintf('\n[%s] Threshold %d / %d : RMSD = %.3f\n', ...
        datestr(now,'HH:MM:SS'), k, nT, rmsdThreshold);

    %% ---- Build adjacency matrix ----
    t_adj = tic;

    RMSD_binary = RMSD <= rmsdThreshold;

    fprintf('    Adjacency built in %.2f s\n', toc(t_adj));

    %% ---- Build graph ----
    t_graph = tic;

    G = graph(RMSD_binary, 'omitselfloops', 'upper');

    fprintf('    Graph construction in %.2f s | Nodes = %d | Edges = %d\n', ...
        toc(t_graph), numnodes(G), numedges(G));

    %% ---- Centrality calculation ----
    t_cent = tic;

    closenessCent = centrality(G, 'closeness');

    fprintf('    Closeness centrality in %.2f s\n', toc(t_cent));

    %% ---- Metrics ----
    nEdges(k) = numedges(G);   % already undirected, no need /2
    okayness(k) = mean(closenessCent) / sqrt(nEdges(k));

    fprintf('    nEdges = %d | mean(closeness) = %.4e | okayness = %.4e\n', ...
        nEdges(k), mean(closenessCent), okayness(k));

    %% ---- ETA estimate ----
    elapsed = toc(t_total);
    avg_per_k = elapsed / k;
    remaining = avg_per_k * (nT - k);

    fprintf('    Elapsed: %.1f s | ETA: %.1f s (~%.1f min)\n', ...
        elapsed, remaining, remaining/60);

end

fprintf('\n------------------------------------------------------------\n');
fprintf('[%s] Finished threshold sweep in %.2f seconds (%.2f minutes)\n', ...
    datestr(now,'HH:MM:SS'), toc(t_total), toc(t_total)/60);

%% ---- Plotting ----

nexttile([2 1])
box on; grid on
set(gcf,'color','w')
set(gca,'LineWidth',2,'FontSize',15)
plot(rmsdThresholdSearchRange, nEdges, 'ko-', 'LineWidth', 2)
xlim([rmsdThresholdSearchRange(1), rmsdThresholdSearchRange(end)])
xlabel('RMSD threshold (\AA)','Interpreter','latex','FontSize',22)
ylabel('# of edges')
legend('# of edges','Location','NorthWest')

nexttile([2 1])
box on; grid on
set(gcf,'color','w')
set(gca,'LineWidth',2,'FontSize',15)
plot(rmsdThresholdSearchRange, okayness, 'o-', 'LineWidth', 2)
xlim([rmsdThresholdSearchRange(1), rmsdThresholdSearchRange(end)])
xlabel('RMSD threshold (\AA)','Interpreter','latex','FontSize',22)
ylabel('Okayness','FontSize',22)

%% ---- Pick best threshold ----
[okayness_max, whereMax] = max(okayness);
rmsd_threshold = rmsdThresholdSearchRange(whereMax);

fprintf('\n[%s] Optimal RMSD threshold = %.3f (okayness = %.4e)\n', ...
    datestr(now,'HH:MM:SS'), rmsd_threshold, okayness_max);

end



