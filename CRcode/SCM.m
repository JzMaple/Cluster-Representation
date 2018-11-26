% Incremental multi-graph matching
% rawMat: initial global matching n(N+1)*n(N+1), n: node size, N: graph
% rawMat contains only the initial matching within each cluster, and it's
% the composition of previous IMGM result and the new graph
% affMat: edge affinity/similarity
% number, the N+1th is new graph
% affScore: affinity score matrix for each pair of graphs, size N*N
% param: structure storing extra parameters
%        param.method: partition method (1: AP, 2: DPP, 3: random)

function [M,idx] = SCM(affScore, rawMat, param, preIdx, srcore, consistency)
    global affinity
    [a,b]=size(affScore);
    if a~=b
        error('The similarity matrix must be square!');
    end
    n=param.n; N=param.N+1;
    massOutlierMode = 0;
    nCluster = ceil(N / 15 - 0.99);
    if nCluster < 1 
        nCluster = 1; 
    end
    inlierMask =  zeros(n,N);
    nDivide = ones([1 N])*n;
    
    Matching = mat2cell(rawMat,nDivide,nDivide);% divide Matching into blocks
    maxScore = max(affScore(:));
    Similarity = arrayfun(@(x) x/maxScore,affScore);
    Similarity = (Similarity*1).^2;
    medianSimilarity = median(Similarity(:));
    p = ones([1 N])*medianSimilarity;
    %%%%%%%%%%%% visualize graphs if specified %%%%%%%%%%%%%%%%
    if param.visualization == 1
        maxSim = max(Similarity(:));
        tmpSimilarity = Similarity + eye(N)*1.05*maxSim;
        tmpSimilarity = ones(N,N)*1.05*maxSim - tmpSimilarity;
        pointMDS = mdscale(tmpSimilarity.^2, 2);
    end
    isKmeans = 1;
    %%%%%%%%%%%%%%%%% generating the topology of hypergraph %%%%%%%%%%%%%%%
    if param.method == isKmeans
        %tic;
        idx = kMeans(Similarity,nCluster,preIdx);
        C_index = unique(idx);
        %toc
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     CAO multi-graph matching on each cluster          %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nCluster = length(C_index);
    ClusterPosition = cell(nCluster,1);
    ClusterSize = zeros(nCluster,1);
    for i = 1:nCluster
        tmpClusterPosition = find(idx == C_index(i));           % find the indices of graphs which belongs to i-th cluster
        ClusterPosition{i} = tmpClusterPosition;
        iClusterSize = length(tmpClusterPosition);              % iClusterSize records the i-th cluster's size 
        ClusterSize(i) = iClusterSize;                          % ClusterSize records the sizes of all clusters
        affinityCrop = cropAffinity(tmpClusterPosition);        % crop the affinity of i-th cluster
        targetCrop = cropTarget(tmpClusterPosition);            % crop the target of i-th cluster
        tmpRawMat = cell2mat(Matching(tmpClusterPosition, tmpClusterPosition));
        affScoreCurrent = affScore(tmpClusterPosition,tmpClusterPosition);
        scrDenom = max(max(affScoreCurrent(1:end,1:end)));
        tmpMatching = CAO_local(tmpRawMat, n, iClusterSize, param.iterMax, scrDenom, affinityCrop, targetCrop, 'exact', 1);      % perform CAO_pc
        iDivide = ones([1 iClusterSize]) * n;
        Matching(tmpClusterPosition, tmpClusterPosition) = mat2cell(tmpMatching,iDivide,iDivide);                               % record the i-th cluster's matching result in Machting
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Choose Representation by finding k graphs with      %%%
    %%%     max consistency in each cluster. And match them     %%%
    %%%     together by multi-graph matching algorithm          %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nCluster > 1
        for i = 1:nCluster
            for j = 1:nCluster
                iCluster = ClusterPosition{i};
                jCluster = ClusterPosition{j};
                r = Representation(iCluster,jCluster,srcore,consistency,N);
                for x = iCluster
                    for y = jCluster
                        Matching{x,y} = Matching{x,r(1)} * Matching{r(1),r(2)} * Matching{r(2),y};
                        Matching{y,x} = Matching{y,r(2)} * Matching{r(2),r(1)} * Matching{r(1),x};
                    end
                end
            end
        end
    end
    
    M = cell2mat(Matching);
end