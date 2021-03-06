% nInlier: number of inliers
% testk: a number for recording the current test number
function affinity = generateRandomAffinity(nInlier,testk,n)
global target
affinity.BiDir = target.config.affinityBiDir;
affinity.edgeAffinityWeight = target.config.edgeAffinityWeight;
affinity.angleAffinityWeight = target.config.angleAffinityWeight;

bGraphMatch = target.config.bGraphMatch;
bUnaryEnable = target.config.bUnaryEnable;
bEdgeEnable = target.config.bEdgeEnable;
graphCnt = target.config.graphCnt;
Sacle_2D = target.config.Sacle_2D;
deform = target.config.deform;
density = target.config.density;
nOutlier = target.config.nOutlier;
complete = target.config.complete;
nodeCnt = nInlier + nOutlier;
affinity.graphCnt = graphCnt;
affinity.nodeCnt = nodeCnt;
affinity.EG = cell(graphCnt,1);
target.pairwiseMask{testk} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);

basePoint = rand(nInlier,1);
if bGraphMatch% direct generate random edges
    baseEdge = tril(rand(nInlier),-1); baseEdge = baseEdge + baseEdge';
else% first randomly generate the coordinates of 2-D points, then derive the edge features by their distance
    basePointSet = randn(2, nInlier);
end

if complete<1
    target.pairwiseMask{testk} = generateRandomCompleteMask(nodeCnt,graphCnt,target.config.complete);
else % fully complete, hence all are ones
    target.pairwiseMask{testk} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);
end

if target.random && ~bGraphMatch
    trans = cell(n,1);
    scale = cell(n,1);
    for i = 1:n
        trans{i} = deform * randn(2,nInlier);
        lambda_1 = 1 + randn()*0.15;
        lambda_2 = 1 + randn()*0.15;
        scale{i} = [lambda_1 0 ; 0 lambda_2];
    end
end

% perm = randperm(graphCnt);
perm = 1:graphCnt;
for gc=1:graphCnt
    i = floor((gc * n - 0.001) / graphCnt) + 1; % generate by i-th core graph
    affinity.nP{gc} = nodeCnt;
    % MPM (Cho et al. CVPR2014) method is more suitable in the random point set matching case i.e. bCreateRandom=0
    % Hence to fully release the capability of MPM, here we also involve bCreateRandom=0
    if bGraphMatch% randomly and directly generate the edge weight of each graph, without setting coordinates of the points
        randomEdge = tril(rand(nodeCnt),-1);
        % generate the Gaussian noise
        N = deform*tril(randn(nodeCnt),-1);N=N+N';
        affinity.edge{gc} = randomEdge+randomEdge';
        affinity.edge{gc}(1:nInlier,1:nInlier) = baseEdge;% + distortAdd;%baseEdge
        % add the noise
        affinity.edge{gc} = affinity.edge{gc} + N;
        affinity.edgeRaw{gc} = affinity.edge{gc};
        % delete the edges below a density threshold to sparsify the random graph
        P = tril(rand(nodeCnt),-1); P = P+P';
        affinity.adj{gc} = logical(P<density);
        affinity.edge{gc}(P>=density) = NaN;
        % generate random point-wise feature
        randomPoint = rand(1,nodeCnt);
        PN = deform*randn(1,nodeCnt);
        affinity.pointFeat{gc} = randomPoint;
        affinity.pointFeat{gc}(1:nInlier) = basePoint{i};
        affinity.pointFeat{gc} = affinity.pointFeat{gc} + PN;% add some noises
    else% point matching, first generate the coordinate of 2D points, then derive their distance as edge weiths 
        if target.random
            pointTmp = scale{i} * basePointSet + trans{i};
            pointSet = [pointTmp + deform*randn(2,nInlier) randn(2,nOutlier)];
        else
            pmax = max(target.data{perm(gc)},[],2);
            pmin = min(target.data{perm(gc)},[],2);
            pointOut = rand(2,nOutlier) .* (pmax - pmin) + pmin;
            pointSet = [target.data{perm(gc)} pointOut];
        end
        pointSet = pointSet';
        % 2nd Order Matrix
        G = zeros(nodeCnt,nodeCnt);
        for r = 1 : length(pointSet)
            for c = r+1 : length(pointSet)
                G(r,c) = sqrt((pointSet(r,1)-pointSet(c,1))^2+(pointSet(r,2)-pointSet(c,2))^2);
            end
        end
        G = G / max(G(:,:));
        G = G + G';
        affinity.edge{gc} = G;
        affinity.edgeRaw{gc} = G;
        
        if target.fc
            P = tril(rand(nodeCnt),-1); P = P+P';
            affinity.adj{gc} = logical(P<density);
            affinity.edge{gc}(P>=density) = NaN;
        else
            tri = delaunay(pointSet(:,1),pointSet(:,2));
            triNum=size(tri,1);
            affinity.adj{gc} = zeros( nodeCnt, nodeCnt);
            for i=1:triNum
                affinity.adj{gc}(tri(i,1),tri(i,2))=1;
                affinity.adj{gc}(tri(i,2),tri(i,1))=1;
                affinity.adj{gc}(tri(i,2),tri(i,3))=1;
                affinity.adj{gc}(tri(i,3),tri(i,2))=1;
                affinity.adj{gc}(tri(i,1),tri(i,3))=1;
                affinity.adj{gc}(tri(i,3),tri(i,1))=1;
            end
            affinity.edge{gc}(affinity.adj{gc}==0) = NaN;
        end
        
        affinity.pointFeat{gc} = rand(1,nodeCnt);
        affinity.pointFeat{gc}(1:nInlier) = randn(nInlier,1);
        affinity.pointFeat{gc} = affinity.pointFeat{gc} + deform*randn(1,nodeCnt);
    end
    
    affinity.nE{gc} = sum(sum(~isnan(affinity.edge{gc})));
    [r,c]=find(~isnan(affinity.edge{gc}));
    affinity.EG{gc}=[r,c]';% size is 2 \times Data.nE{1} (i.e. number of edges), the fist row is one ending point, the second row is the other
    affinity.edgeFeat{gc} = affinity.edge{gc}(~isnan(affinity.edge{gc}))';% edgeFeat is of size 1\times number of edges
    
    % incidence matrix
    affinity.G{gc} = zeros(nodeCnt,affinity.nE{gc});
    for c = 1 : affinity.nE{gc}
        affinity.G{gc}(affinity.EG{gc}(:, c), c) = 1;
    end
    % augumented adjacency
    affinity.H{gc} = [affinity.G{gc}, eye(nodeCnt)];
end
for xview = 1:graphCnt
    if affinity.BiDir
        yviewSet = [1:xview-1,xview+1:graphCnt];
    else
        yviewSet = xview+1:graphCnt;
    end
    for yview = yviewSet
        if bUnaryEnable
            affinity.KP{xview,yview} = exp(-conDst(affinity.pointFeat{xview}, affinity.pointFeat{yview},0) / Sacle_2D);
        else
            affinity.KP{xview,yview} = zeros(nodeCnt,nodeCnt);
        end
        if bEdgeEnable
            affinity.KQ{xview,yview} = exp(-conDst(affinity.edgeFeat{xview}, affinity.edgeFeat{yview},0) / Sacle_2D);%但边和边之间有affinity Kq，用指数形式比较平滑便于优化,形成指数化后的边边矩阵n1*n2
        else
            affinity.KQ{xview,yview} = zeros(nodeCnt^2,nodeCnt^2);
        end
        affinity.K{xview,yview} = conKnlGphKU(affinity.KP{xview,yview}, affinity.KQ{xview,yview}, affinity.EG{xview},affinity.EG{yview});%EG是边的端点的索引，2*n1,2*n
        affinity.K{xview,yview} = full(affinity.K{xview,yview});
    end
end