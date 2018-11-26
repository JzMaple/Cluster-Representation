clear all; clc; close all;
isShowFig = 1; isShowTri = 1; isShowResult = 1;
isCombineTri = 0; % if we combine triangulation of two partitions
delta = 150; delta2 = 200000;
nTestDensity = 10;
imgDir = 'cutmatch_image/';
imgFiles = dir([imgDir '*.jpg']);
testDensity = [0:0.02:0.04];
isLength = 0; % whether to use the distance and connectivity for cuts

for fptr = 1:length(imgFiles)
    display([num2str(fptr)]);
    img = imread([imgDir imgFiles(fptr).name]);
    img = img(144:455,51:677,:);
    imgGray = rgb2gray(img);
    [m, n]=size(imgGray);
    tmpName = imgFiles(fptr).name;
    imgLetter = tmpName(1:(find(tmpName=='.')-1));
    fp = fopen([imgDir imgLetter '.txt'],'r');
    tmpLine = fgets(fp);
    nNode = str2num(tmpLine);
    partition1 = zeros(2,nNode);
    partition2 = zeros(2,nNode);
    for j = 1:nNode
        tmpLine = fgets(fp);
        a = sscanf(tmpLine,'%d %d %d %d \n');
        %partition1(1,j) = a(1); partition1(2,j) = a(2);
        %partition2(1,j) = a(3); partition2(2,j) = a(4);
        partition1(1,j) = a(1)-51; partition1(2,j) = a(2)-144;
        partition2(1,j) = a(3)-51; partition2(2,j) = a(4)-144;
    end
    fclose(fp);
    %partition1(1,:) = m -partition1(1,:);
    partition = [partition1, partition2];
    partition =  partition';
    if isShowFig
        figure; imshow(img); hold on;
        scatter(partition1(1,:),partition1(2,:),40,'r','MarkerFaceColor','r');
        scatter(partition2(1,:),partition2(2,:),40,'b','MarkerFaceColor','b');
        hold off;
    end
    %%%%% delaunay on partition 1 and 2 %%%%%
    tri1 = delaunayTriangulation(partition1');
    tri2 = delaunayTriangulation(partition2');
    triEdge1 = edges(tri1);
    triEdge2 = edges(tri2);
    %%%%% delaunay on whole graph %%%%%
    tri = delaunayTriangulation([partition1,partition2]');
    triEdge = edges(tri);
    GTmatching = [zeros(nNode,nNode),eye(nNode);eye(nNode),zeros(nNode,nNode)];
    if isShowTri
        figure; imshow(img); hold on; box off; title('Ground-truth partitions and matchings','FontSize',18);
        triplot(tri.ConnectivityList,[partition1(1,:),partition2(1,:)],[partition1(2,:),partition2(2,:)]);
        scatter(partition1(1,:),partition1(2,:),40,'r','MarkerFaceColor','r');
        scatter(partition2(1,:),partition2(2,:),40,'b','MarkerFaceColor','b');
        for i = 1:2*nNode-1
            for j = i:2*nNode
                if GTmatching(i,j)==1
                    line([partition(i,1),partition(j,1)],[partition(i,2),partition(j,2)],'LineStyle','--','Color','g','LineWidth',2);
                end
            end
        end
    end
    %%%%% generate mask for baseline graph %%%%%
    [X, Y]=meshgrid(1:2*nNode,1:2*nNode);
    LL = [X(:),Y(:)];
    G = partition(LL(:,1),:)-partition(LL(:,2),:);
    G = sqrt(G(:,1).^2+G(:,2).^2);
    G = reshape(G, [2*nNode 2*nNode]); % distance between each pair of nodes
    mask = zeros(2*nNode,2*nNode);
    for k = 1:size(triEdge,1)
        mask(triEdge(k,1),triEdge(k,2))=1;
    end
    if isCombineTri
        for k = 1:size(triEdge1,1)
            mask(triEdge1(k,1),triEdge1(k,2))=1;
            mask(triEdge1(k,1)+nNode,triEdge1(k,2)+nNode)=1;
        end
        for k = 1:size(triEdge2,1)
            mask(triEdge2(k,1)+nNode,triEdge2(k,2)+nNode)=1;
            mask(triEdge2(k,1),triEdge2(k,2))=1;
        end
    end
    mask = mask + mask'; % edge mask
    
    nBaseline = length(find(mask == 1))/2;
    nFullCon = 2*nNode*(2*nNode-1)/2;
    nExtraNode = round((nFullCon - nBaseline)*testDensity);
    
    densityList = [];
    for i=1:2*nNode-1
        for j = i+1:2*nNode
            if mask(i,j) == 0
                densityList = [densityList,[i,j]'];
            end
        end
    end
    
    %%%%%% generate sift des for each node %%%%%%
    nodeFeature = zeros(128,2*nNode);
    for i = 1:2*nNode
        fc = [partition(i,:)'; 5; 0];
        [f,des]=vl_sift(single(imgGray),'frames',fc,'orientations');
        nodeFeature(:,i) = double(des(:,1));
        if isLength
            
        end
    end

    %%%%% loop for varying edge density %%%%%
    for iDensityPtr = 1:length(testDensity)
        %%%%% multiple exp for each density %%%%%
        for iTest = 1:nTestDensity
            %%%%% generating edge mask for varying density %%%%%
            if testDensity(iDensityPtr) == 0
                maskCur = mask;
                if iTest>1
                    acc_matching{fptr,iDensityPtr,iTest} = acc_matching{fptr,iDensityPtr,1};
                    acc_cut{fptr,iDensityPtr,iTest} = acc_cut{fptr,iDensityPtr,1};
                    acc_matching_rrwm{fptr,iDensityPtr,iTest} = acc_matching_rrwm{fptr,iDensityPtr,1};
                    acc_matching_IBGP{fptr,iDensityPtr,iTest} = acc_matching_IBGP{fptr,iDensityPtr,1};
                    acc_matching_IPFP{fptr,iDensityPtr,iTest} = acc_matching_IPFP{fptr,iDensityPtr,1};
                    acc_cut_raw{fptr,iDensityPtr,iTest} = acc_cut_raw{fptr,iDensityPtr,1};
                    continue;
                end
            elseif testDensity(iDensityPtr) == 1
                maskCur = ones(2*nNode,2*nNode) - eye(2*nNode);
            else
                maskTmp = zeros(2*nNode,2*nNode);
                randList = randperm(size(densityList,2));
                nSelect = round(testDensity(iDensityPtr)*size(densityList,2));
                for i = 1:nSelect
                    maskTmp(densityList(1,randList(i)),densityList(2,randList(i)))=1;
                end
                maskTmp = maskTmp + maskTmp';
                maskCur = mask + maskTmp;
            end
            
            [X, Y]=meshgrid(1:(2*nNode)^2,1:(2*nNode)^2);
            X = X(:);  Y = Y(:);
            G = G(:);
            
            aff = exp(-(G(X) - G(Y)).^2/delta);
            affinityA = reshape(aff, [((2*nNode)^2),((2*nNode)^2)]); % initialize affinity matrix
            
            maskE = zeros([((2*nNode)^2),((2*nNode)^2)]);
            for i = 1:2*nNode
                for j = 1:2*nNode;
                    for k = 1:2*nNode
                        for l = 1:2*nNode
                            if maskCur(i,j)==1 && maskCur(k,l)==1 && (i~=k || j~=l)
                                maskE(sub2ind([2*nNode, 2*nNode],i,k),sub2ind([2*nNode, 2*nNode],j,l))=1;
                            end
                        end
                    end
                end
            end
            affinityA = affinityA.*maskE; % sparse affinity A according to graph connectivity
            % continue;
            %%%%%% weight for cut %%%%%
            featureWeight = zeros(2*nNode,2*nNode);
            for i = 1:2*nNode
                for j=1:2*nNode
                    if maskCur(i,j)==1
                        featureWeight(i,j)=exp(-(norm(nodeFeature(:,i)-nodeFeature(:,j)))^2/delta2); % kernel function: similarity
                        if isLength
                            featureWeight(i,j) = featureWeight(i,j) + exp(-(norm(partition(i,:)-partition(j,:)))^2/500);
                            %featureWeight(i,j) =  exp(-norm(partition(i,:)-partition(j,:))^2/500);
                        end
                        % affinity.featureWeight(i,j)=2*(1/(1+exp(-abs(0.01*norm(affinity.nodeFeature(:,i)-affinity.nodeFeature(:,j),1))))-0.5); % sigmoid function: dissimilarity
                    end
                end
            end
            
            %%%%% try raw graph cut %%%%%
            D = diag(sum(featureWeight,1));
            Lap = D - featureWeight;
            [eigVec,eigVal] = eigs(Lap,2,'SA');
            GClabel = sign(eigVec(:,2));
            
            if isShowResult
                figure; imshow(img); hold on; axis off; box off; title('Initial graph cut result','FontSize',18);
                triplot(tri.ConnectivityList,[partition1(1,:),partition2(1,:)],[partition1(2,:),partition2(2,:)]);
                scatter(partition(find(GClabel==1),1),partition(find(GClabel==1),2),40,'r','MarkerFaceColor','c');
                scatter(partition(find(GClabel==-1),1),partition(find(GClabel==-1),2),40,'r','MarkerFaceColor','y');
                hold off;
            end
            
            %%%%% try raw RRWM %%%%%
            algpar.iterMax1 = 300;
            X_rrwm_init = 1/(2*nNode)*ones(2*nNode,2*nNode);
            X_rrwm = RRWM(affinityA,2*nNode, 2*nNode,algpar);
            maxWeight2 = max(X_rrwm(:));
            X_rrwm = reshape(X_rrwm,[2*nNode,2*nNode]);
            [MatchingRRWM cost] = hungarian(maxWeight2*ones(2*nNode,2*nNode)-X_rrwm);
            
            %%%%% try raw IPFP %%%%%
            X_ipfp = IPFP(affinityA,2*nNode,2*nNode);
            maxWeightIPFP = max(X_ipfp(:));
            X_ipfp = reshape(X_ipfp,[2*nNode,2*nNode]);
            [MatchingIPFP cost] = hungarian(maxWeightIPFP*ones(2*nNode,2*nNode)-X_ipfp);
            
            %%%%% try IBGP %%%%%
            X_init = ones([2*nNode 2*nNode])/(2*nNode);
            X_init = BregmanBiStochZero(X_init,50);
            X_init = X_init(:);
            
            %tic;
            X_next = IBGP(X_init, affinityA, zeros(2*nNode,1), 0, 0.002, 200, 100, 1);
            %toc;
            X_out = reshape(X_next, [2*nNode, 2*nNode]);
            
            maxWeight = max(X_out(:));
            [matching_initial cost] = hungarian(maxWeight*ones(2*nNode,2*nNode)-X_out);
            
            if isShowResult
                figure; imshow(img); hold on; axis off; box off; title('Matching result with RRWM','FontSize',18);
                triplot(tri.ConnectivityList,[partition1(1,:),partition2(1,:)],[partition1(2,:),partition2(2,:)]);
                scatter(partition1(1,:),partition1(2,:),40,'r','MarkerFaceColor','r');
                scatter(partition2(1,:),partition2(2,:),40,'b','MarkerFaceColor','b');
                for i = 1:2*nNode-1
                    for j = i:2*nNode
                        if MatchingRRWM(i,j)==1
                            line([partition(i,1),partition(j,1)],[partition(i,2),partition(j,2)],'LineStyle','--','Color','g','LineWidth',2);
                        end
                    end
                end
                hold off;
            end
            
            if isShowResult
                figure; imshow(img); hold on; axis off; box off; title('Matching result with IBGP','FontSize',18);
                triplot(tri.ConnectivityList,[partition1(1,:),partition2(1,:)],[partition1(2,:),partition2(2,:)]);
                scatter(partition1(1,:),partition1(2,:),40,'r','MarkerFaceColor','r');
                scatter(partition2(1,:),partition2(2,:),40,'b','MarkerFaceColor','b');
                for i = 1:2*nNode-1
                    for j = i:2*nNode
                        if matching_initial(i,j)==1
                            line([partition(i,1),partition(j,1)],[partition(i,2),partition(j,2)],'LineStyle','--','Color','g','LineWidth',2);
                        end
                    end
                end
                hold off;
            end
            
            %%%%% param of CutMatch -200-50 %%%%%
            param.iterCnt = 5;
            param.lambda_1 = 300;
            param.lambda_2 = 60; 
            param.typeL = 1;
            param.X_initial = 1/(2*nNode)*ones(2*nNode,2*nNode);
            
            %%%%% apply CutMatch %%%%%
            %tic;
            [Matching, Cut]=CutMatching(affinityA, featureWeight,param);
            %toc;
            tmpLabel = sign(Cut);
            maxWeight = max(Matching(:));
            [matching cost]=hungarian(maxWeight*ones(2*nNode,2*nNode)-Matching);
            
            if isShowResult
                figure; imshow(img); hold on; axis off; box off; title('Matching result with CutMatch','FontSize',18);
                triplot(tri.ConnectivityList,[partition1(1,:),partition2(1,:)],[partition1(2,:),partition2(2,:)]);
                scatter(partition1(1,:),partition1(2,:),40,'r','MarkerFaceColor','r');
                scatter(partition2(1,:),partition2(2,:),40,'b','MarkerFaceColor','b');
                for i = 1:2*nNode-1
                    for j = i:2*nNode
                        if matching(i,j)==1
                            line([partition(i,1),partition(j,1)],[partition(i,2),partition(j,2)],'LineStyle','--','Color','g','LineWidth',2);
                        end
                    end
                end
                hold off;
            end
            
            if isShowResult
                figure; imshow(img); axis off; box off; title('Cut result with CutMatch','FontSize',18);
                hold on; triplot(tri.ConnectivityList,[partition1(1,:),partition2(1,:)],[partition1(2,:),partition2(2,:)]);
                scatter(partition(find(tmpLabel==1),1),partition(find(tmpLabel==1),2),40,'r','MarkerFaceColor','w');
                scatter(partition(find(tmpLabel==-1),1),partition(find(tmpLabel==-1),2),40,'r','MarkerFaceColor','k'); hold off;
            end
            
            a = 1;
            %%%%% calculate accuracy now %%%%%
            GTmatching = [zeros(nNode,nNode),eye(nNode);eye(nNode),zeros(nNode,nNode)];
            GTcut = [ones(nNode,1);-ones(nNode,1)];
            diff_matching = GTmatching.*matching;
            diff_matching_rrwm = GTmatching.*MatchingRRWM;
            diff_matching_IPFP = GTmatching.*MatchingIPFP;
            diff_matching_IBGP = GTmatching.*matching_initial;
            diff_cut_num = max([sum(tmpLabel==GTcut),sum(tmpLabel~=GTcut)]);
            diff_cut_numraw = max([sum(GClabel==GTcut),sum(GClabel~=GTcut)]);
            
            acc_matching{fptr,iDensityPtr,iTest} = length(find(diff_matching==1))/(2*nNode);
            acc_cut{fptr,iDensityPtr,iTest} = diff_cut_num/(2*nNode);
            acc_matching_rrwm{fptr,iDensityPtr,iTest} = length(find(diff_matching_rrwm==1))/(2*nNode);
            acc_matching_IBGP{fptr,iDensityPtr,iTest} = length(find(diff_matching_IBGP==1))/(2*nNode);
            acc_matching_IPFP{fptr,iDensityPtr,iTest} = length(find(diff_matching_IPFP==1))/(2*nNode);
            acc_cut_raw{fptr,iDensityPtr,iTest} = diff_cut_numraw/(2*nNode);
            a = 1;
        end
    end
end


for j = 1:length(testDensity)
    acc_ave_cutmatch_gm{j}=0;
    acc_ave_cutmatch_gc{j}=0;
    acc_ave_rrwm_gm{j}=0;
    acc_ave_ipfp_gm{j}=0;
    acc_ave_rawcut_gc{j}=0;
    acc_ave_ibgp_gm{j}=0;
    for i = 1:length(imgFiles)
        for k = 1:nTestDensity
            acc_ave_cutmatch_gm{j} = acc_ave_cutmatch_gm{j} + acc_matching{i,j,k};
            acc_ave_cutmatch_gc{j} = acc_ave_cutmatch_gc{j} + acc_cut{i,j,k};
            acc_ave_rrwm_gm{j} = acc_ave_rrwm_gm{j} + acc_matching_rrwm{i,j,k};
            acc_ave_rawcut_gc{j} = acc_ave_rawcut_gc{j} + acc_cut_raw{i,j,k};
            acc_ave_ibgp_gm{j} = acc_ave_ibgp_gm{j} + acc_matching_IBGP{i,j,k};
        end
    end
    acc_ave_cutmatch_gm{j} = acc_ave_cutmatch_gm{j}/(length(imgFiles)*nTestDensity);
    acc_ave_cutmatch_gc{j} = acc_ave_cutmatch_gc{j}/(length(imgFiles)*nTestDensity);
    acc_ave_rrwm_gm{j} = acc_ave_rrwm_gm{j}/(length(imgFiles)*nTestDensity);
    acc_ave_rawcut_gc{j} = acc_ave_rawcut_gc{j}/(length(imgFiles)*nTestDensity);
    acc_ave_ibgp_gm{j} = acc_ave_ibgp_gm{j}/(length(imgFiles)*nTestDensity);
end

figure; hold on; title('Matching performance on real-world images');
p1 = plot(testDensity,cell2mat(acc_ave_cutmatch_gm),'LineWidth',2,'LineStyle','-','Marker','o','MarkerSize',6,'Color','b');
p2 = plot(testDensity,cell2mat(acc_ave_ibgp_gm),'LineWidth',2,'LineStyle','--','Marker','*','MarkerSize',6,'Color','g');
p3 = plot(testDensity,cell2mat(acc_ave_rrwm_gm),'LineWidth',2,'LineStyle','-','Marker','p','MarkerSize',6,'Color','r');
p4 = plot(testDensity,cell2mat(acc_ave_ipfp_gm),'LineWidth',2,'LineStyle','--','Marker','+','MarkerSize',6,'Color','w');
legend([p1,p2,p3],'CutMatch','IBGP','RRWM');
ylabel('\fontname{times new roman}Accuracy','FontSize',15);
xlabel('\fontname{times new roman}Edge density','FontSize',15);
hold off;

figure; hold on; title('Cuts performance on real-world images');
p1 = plot(testDensity,cell2mat(acc_ave_cutmatch_gc),'LineWidth',2,'LineStyle','-','Marker','o','MarkerSize',6,'Color','b');
p2 = plot(testDensity,cell2mat(acc_ave_rawcut_gc),'LineWidth',2,'LineStyle','--','Marker','*','MarkerSize',6,'Color','g');
legend([p1,p2],'CutMatch','Raw Cuts');
ylabel('\fontname{times new roman}Accuracy','FontSize',15);
xlabel('\fontname{times new roman}Edge density','FontSize',15);
hold off;
