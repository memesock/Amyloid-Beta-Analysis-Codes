function C = determineNKmeans(U)

%% Assign structures within an ensemble to unique classes based on nearest-neighbor RMSD minimization 
%   
%  Determine the optimal # of K-means clusters in a graph representation of
%  a structural ensemble, by maximizing the average silhouette across
%  numerous trial runs. 
%

kRange = 3:1:6;

for N=5
    for k = kRange
        C{N,k} = kmeans(U, k, 'start', 'cluster') - 1;
        s{N,k} = silhouette(U,C{N,k});
    end
end

sAvg{N,k} = mean(s);

end



