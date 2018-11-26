function consistency = cal_cluster_consistency(cluster,rawMat)
    global target
    graphCnt = length(cluster);
    constList = zeros(graphCnt,1);
    nodeCnt = target.config.nodeCnt;
    for ref = 1: graphCnt
        viewk = 1;
        err = zeros((graphCnt+1)*(graphCnt-2)/2,1);
        rscope = (ref-1)*nodeCnt+1:ref*nodeCnt;
        for i = 1:graphCnt
            iscope = (cluster(i)-1)*nodeCnt+1:cluster(i)*nodeCnt;
            for j = i+1:graphCnt
                jscope = (cluster(j)-1)*nodeCnt+1:cluster(j)*nodeCnt;
                Xij = rawMat(iscope,jscope);
                Xirj = X(iscope,rscope)*X(rscope,jscope);
                err(viewk) = sum(sum(abs(Xirj-Xij)))/2/nodeCnt;
                viewk = viewk + 1;
            end
        end
        constList(ref) = 1 - mean(err);
    end
    consistency = mean(constList);
end