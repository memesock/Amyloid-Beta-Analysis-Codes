function varMetricS = plotWireframeMeanVar(S_allX,S_allY,S_allZ)
%Average wireframe positions at each index and plot class averaged
%conformer models, showing spatial variance as a variability metric
%   Prior to calling this function, ensure that the right axes are
%   specified onto which to plot

% Get mean positions at each index for averaged wireframe 
meanS = [mean(S_allX,2), mean(S_allY,2), mean(S_allZ,2)];

% Colour point positions by spatial variance across the set of structures (kind of like a crystallographic B-factor map)
varS = [var(S_allX,0,2), var(S_allY,0,2), var(S_allZ,0,2)];
varMetricS = sqrt(varS(:,1).^2 + varS(:,2).^2 + varS(:,3).^2);
scatter3(meanS(:,1),meanS(:,2),meanS(:,3),300,varMetricS,'filled','MarkerEdgeColor','K','LineWidth',1)
format long g

%a.FontSize = 15;
plot3(meanS(:,1),meanS(:,2),meanS(:,3),'K-','LineWidth',1.5)
plot3(meanS(1,1),meanS(1,2),meanS(1,3),'rs','MarkerFaceColor','red','MarkerSize',20) % show 1st atom
plot3(meanS(end,1),meanS(end,2),meanS(end,3),'r^','MarkerFaceColor','red','MarkerSize',10) % show last atom
box on
axis off

end