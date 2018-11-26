function idx = kMeans(Similarity,nCluster,preIdx)
    n = size(Similarity,1);
    id = [preIdx 1];
    opt = cal_cluster_value(id,Similarity,nCluster);
    flag = 1;
    while flag
        flag = 0;
        for grh = 1 : n % for each graph, find the best cluster it belongs to
            optCluster = id(grh);
            for cl = 1 : nCluster
                id(grh) = cl;
                tmp = cal_cluster_value(id,Similarity,nCluster);
                if tmp > opt
                    opt = tmp;
                    optCluster = cl; 
                end
            end
            id(grh) = optCluster;
        end
    end
    idx = id;
end