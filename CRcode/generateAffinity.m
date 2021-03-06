function [affinity]= generateAffinity()
global target GT
bUnaryEnable = target.config.bUnaryEnable;
bEdgeEnable = target.config.bEdgeEnable;
affinity.BiDir = target.config.affinityBiDir;
% affinity.Factorize = target.config.affinityFGM;
% adjFlag = target.config.adjFlag;
affinity.edgeAffinityWeight = target.config.edgeAffinityWeight;
affinity.angleAffinityWeight = target.config.angleAffinityWeight;
Sacle_2D = target.config.Sacle_2D;
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;
affinity.graphCnt = graphCnt;
affinity.nodeCnt = nodeCnt;
adjlen = zeros(graphCnt,1);
Data = cell(graphCnt,1);
for viewk = 1:graphCnt
    Data{viewk}.nP = size(target.data{viewk}.point,1);%点的数目
    Data{viewk}.edge = zeros(nodeCnt,nodeCnt);%边
    Data{viewk}.point = target.data{viewk}.point;%点的坐标
    Data{viewk}.angle = zeros(Data{viewk}.nP,Data{viewk}.nP);%角度
%     %考虑一阶sift特征
%     for s=1:nodeCnt
%         Data{viewk}.edge(s,s) = data{viewk}.feat
%     end
    %计算每条边的长度edge和角度angle    
    for r = 1:nodeCnt
        for c = r+1:nodeCnt
            Data{viewk}.edge(r,c) = sqrt((target.data{viewk}.point(r,1)-target.data{viewk}.point(c,1))^2+(target.data{viewk}.point(r,2)-target.data{viewk}.point(c,2))^2);
            % atan(X) is in the range [-pi/2,pi/2]. 现在再转成[-90,90]
            Data{viewk}.angle(r,c) = 180/pi*atan((target.data{viewk}.point(r,2)-target.data{viewk}.point(c,2))/(target.data{viewk}.point(r,1)-target.data{viewk}.point(c,1)));
        end
    end
    Data{viewk}.edge = Data{viewk}.edge/max(Data{viewk}.edge(:));
    Data{viewk}.edge = Data{viewk}.edge + Data{viewk}.edge';
    % 再转成[-1,1]
    Data{viewk}.angle = Data{viewk}.angle/90;
    Data{viewk}.angle = Data{viewk}.angle + Data{viewk}.angle';
            
    Data{viewk}.edgeRaw = Data{viewk}.edge;
    Data{viewk}.angleRaw = Data{viewk}.angle;
    
    if adjFlag>1%1是full connect
        tri = delaunay(Data{viewk}.point(:,1),Data{viewk}.point(:,2));
        triNum=size(tri,1);
        Data{viewk}.adjMatrix = zeros( Data{viewk}.nP, Data{viewk}.nP);
        for i=1:triNum
            Data{viewk}.adjMatrix(tri(i,1),tri(i,2))=1;
            Data{viewk}.adjMatrix(tri(i,2),tri(i,1))=1;
            Data{viewk}.adjMatrix(tri(i,2),tri(i,3))=1;
            Data{viewk}.adjMatrix(tri(i,3),tri(i,2))=1;
            Data{viewk}.adjMatrix(tri(i,1),tri(i,3))=1;
            Data{viewk}.adjMatrix(tri(i,3),tri(i,1))=1;
        end
    else%1是full connect
        Data{viewk}.adjMatrix = ones(Data{viewk}.nP, Data{viewk}.nP);
    end
    adjlen(viewk) = sum(Data{viewk}.adjMatrix(:));
end
switch adjFlag
    case {1,5}% 1:full connected,5:各管各
        for viewk=1:graphCnt
            Data{viewk}.adjMatrix = logical(Data{viewk}.adjMatrix);
            Data{viewk}.nE = sum(Data{viewk}.adjMatrix(:));
        end
    case {2,3,4}% 2:并集,3:最大集,4:交集,当顺序被打乱的时候，单纯取最多边数不work 需要先重对齐 再取最大边数
        maxid = find(adjlen==max(adjlen));reference=maxid(1);
        refAdj = Data{reference}.adjMatrix;
        
        if adjFlag==2||adjFlag==4%2是并集，4是交集，需要进一步充实或者稀疏reference adj
            for viewk=1:graphCnt
                if viewk~=reference
                    %更新adjMatrix: reference (r) vs. viewk (v)
                    permGT = GT((reference-1)*nodeCnt+1:reference*nodeCnt,(viewk-1)*nodeCnt+1:viewk*nodeCnt);
                    permgt = mat2perm(permGT);%permgt(r)=v
                    if adjFlag==2%并集 做加法
                        for i=1:nodeCnt
                            for j=i+1:nodeCnt
                                refAdj(i,j) = refAdj(i,j)+Data{viewk}.adjMatrix(permgt(i),permgt(j));
                            end
                        end
                    elseif adjFlag==4%交集 做乘法
                        for i=1:nodeCnt
                            for j=i+1:nodeCnt
                                refAdj(i,j) = refAdj(i,j)*Data{viewk}.adjMatrix(permgt(i),permgt(j));
                            end
                        end
                    end
                    refAdj = refAdj + refAdj';
                end
            end
        end
        refAdj = logical(refAdj);
        %下面用reference adj来更新其他adj
        for viewk=1:graphCnt
            if viewk~=reference
                %更新adjMatrix: viewk (v) vs. reference (r)
                permGT = GT((viewk-1)*nodeCnt+1:viewk*nodeCnt,(reference-1)*nodeCnt+1:reference*nodeCnt);
                permgt = mat2perm(permGT);%permgt(v)=r
                adj = zeros(nodeCnt,nodeCnt);
                for i=1:nodeCnt
                    for j=i+1:nodeCnt
                        adj(i,j) = refAdj(permgt(i),permgt(j));
                    end
                end
                Data{viewk}.adjMatrix = adj + adj';
            end
        end
end

affinity.EG = cell(graphCnt,1);

for gc=1:graphCnt
    Data{gc}.adjMatrix = logical(Data{gc}.adjMatrix);
    Data{gc}.edge(~Data{gc}.adjMatrix) = NaN;
    Data{gc}.angle(~Data{gc}.adjMatrix) = NaN;
    Data{gc}.nE = sum(sum(Data{gc}.adjMatrix));
            
    [r,c]=find(~isnan(Data{gc}.edge));
    affinity.EG{gc}=[r,c]';%2*Data.nE{1} 即2*边的数目 第一行是起点 第二行是终点
    Data{gc}.edgeFeat = Data{gc}.edge(~isnan(Data{gc}.edge))';%edgeFeat是一个1*边数的矩阵,edge是triangle的mask
    Data{gc}.angleFeat = Data{gc}.angle(~isnan(Data{gc}.angle))';
    if bUnaryEnable
        Data{gc}.pointFeat = target.data{gc}.feat';%一列是一个点的feature @todo
    end
    affinity.nP{gc} = Data{gc}.nP;
    % incidence matrix
    affinity.G{gc} = zeros(Data{gc}.nP, Data{gc}.nE);
    for c = 1 : Data{gc}.nE
        affinity.G{gc}(affinity.EG{gc}(:, c), c) = 1;
    end
    % augumented adjacency
    affinity.H{gc} = [affinity.G{gc}, eye(Data{gc}.nP)];
    affinity.edge{gc} = Data{gc}.edge;
    affinity.edgeRaw{gc} = Data{gc}.edgeRaw;%这个可以直接用来做admm
    affinity.angleRaw{gc} = Data{gc}.angleRaw;
    affinity.adj{gc} = Data{gc}.adjMatrix;
    
end
%下面计算affinity矩阵
for xview = 1:graphCnt
%     for yview = xview+1:graphCnt
    if affinity.BiDir
        yviewSet = [1:xview-1,xview+1:graphCnt];
    else
        yviewSet = xview+1:graphCnt;
    end
    for yview = yviewSet
        % % % 考虑一阶特征相似度     
        if bUnaryEnable%@todo
            featAffinity = conDst(Data{xview}.pointFeat, Data{yview}.pointFeat,0)/10000/128;
            affinity.KP{xview,yview} = exp(-featAffinity/ Sacle_2D);%怎么设置权重？？
        else
            affinity.KP{xview,yview} = zeros(Data{xview}.nP, Data{yview}.nP);
        end
        % % %
        dq = zeros(length(Data{xview}.edgeFeat),length(Data{yview}.edgeFeat));
        if bEdgeEnable
            if isfield(Data{xview},'edgeFeat') && affinity.edgeAffinityWeight>0
                dq = dq + affinity.edgeAffinityWeight*conDst(Data{xview}.edgeFeat, Data{yview}.edgeFeat,0);%计算两个图里边边之间的距离平方，形成边的距离矩阵n1*n2
            end
            if isfield(Data{xview},'angleFeat') && affinity.angleAffinityWeight>0
                dq = dq + affinity.angleAffinityWeight*conDst(Data{xview}.angleFeat, Data{yview}.angleFeat,1);
            end
    %         affinity.DQ{xview,yview} = dq;
            affinity.KQ{xview,yview} = exp(-dq / Sacle_2D);%但边和边之间有affinity Kq，用指数形式比较平滑便于优化,形成指数化后的边边矩阵n1*n2
    %         affinity.Ct{xview,yview} = ones(size(affinity.KP{xview,yview}));
        else
            affinity.KQ{xview,yview} = dq;
        end
        affinity.K{xview,yview} = conKnlGphKU(affinity.KP{xview,yview}, affinity.KQ{xview,yview}, affinity.EG{xview},affinity.EG{yview});%EG是边的端点的索引，2*n1,2*n2
    end
end

%下面在affinity上增加一些factorized的信息
% 即factMultiGraph函数的内容
if affinity.Factorize
    for i=1:graphCnt
        affinity.HH{i} = affinity.H{i}' * affinity.H{i};
        affinity.IndG{i} = mat2ind(affinity.G{i});
        affinity.IndGT{i} = mat2ind(affinity.G{i}');
        affinity.IndH{i} = mat2ind(affinity.H{i});
        affinity.IndHT{i} = mat2ind(affinity.H{i}');
    end
    
    affinity.L = cell(graphCnt,graphCnt);
    affinity.UU = cell(graphCnt,graphCnt);
    affinity.VV = cell(graphCnt,graphCnt);
    affinity.UUHH = cell(graphCnt,graphCnt);
    affinity.VVHH = cell(graphCnt,graphCnt);
    affinity.HUH = cell(graphCnt,graphCnt);
    affinity.HVH = cell(graphCnt,graphCnt);
    affinity.GKGP = cell(graphCnt,graphCnt);
    
    
    
    for i=1:graphCnt
        if affinity.BiDir
            jset = [1:i-1,i+1:graphCnt];
        else
            jset = i+1:graphCnt;
        end
        for j=jset
            % L
%             fprintf('i=%d,j=%d\n',i,j);
            affinity.L{i,j} = [affinity.KQ{i,j}, -affinity.KQ{i,j} * affinity.G{j}'; 
                -affinity.G{i} * affinity.KQ{i,j}, affinity.G{i} * affinity.KQ{i,j} * affinity.G{j}' + affinity.KP{i,j}];
            % SVD
            [U, S, V] = svd(affinity.L{i,j});
            s = diag(S);
            idx = 1 : length(s);
            %idx = find(s > 0);
    %         pr('svd: %d/%d', length(idx), length(s));
    %         k = length(idx);
            U = U(:, idx);
            V = V(:, idx);
            s = s(idx);

            U = multDiag('col', U, real(sqrt(s)));
            V = multDiag('col', V, real(sqrt(s)));

            % the following decomposition works very badly
            % U = eye(size(L, 1));
            % V = L;

            % auxiliary variables that will be frequently used in the optimization
            affinity.UU{i,j} = U * U';
            affinity.VV{i,j} = V * V';

            affinity.UUHH{i,j} = affinity.UU{i,j} .* affinity.HH{i};
            affinity.VVHH{i,j} = affinity.VV{i,j} .* affinity.HH{j};
            affinity.HUH{i,j} = affinity.H{i} * affinity.UUHH{i,j} * affinity.H{i}';
            affinity.HVH{i,j} = affinity.H{j} * affinity.VVHH{i,j} * affinity.H{j}';
            affinity.GKGP{i,j} = -affinity.G{i} * affinity.KQ{i,j} * affinity.G{j}' + affinity.KP{i,j};
        end
    end
end

