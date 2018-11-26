function v = cal_cluster_value(idx,Similarity,nCluster)
    value = 0;
    for cl = 1 : nCluster
        cluster = find(idx == cl);
        if isempty(cluster) continue; end
        value = value + mean(mean(Similarity(cluster,cluster))); 
    end
    v = value;
end