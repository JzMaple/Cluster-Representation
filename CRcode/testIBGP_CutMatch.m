clear all; clc; close all; warning off;
testMode = 1;
halfNodeCnt = 20;
m = rand([10 1]);
m = m/sum(m);
nTest = 3;
isDisplay = 0;

param.deform = 0.1; % random noise level of point position
param.gap = -0.1; % overlapping level
param.saveGraph = 1; param.delta = 0.5; param.isTri = 1; % if delaunay triangulation
param.isNodeFeature = 1;
param.featureDeform = 0.1; % level to add noise to node feature
param.isFeatureWeight = 1;
param.iterCnt = 5;
param.lambda_1 = 200; lambda_1 = 1;
param.lambda_2 = 50; lambda_2 = 0.1;
param.typeL = 1;
param.X_initial = 1/(2*halfNodeCnt)*ones(2*halfNodeCnt,2*halfNodeCnt);
% param.variationProb = 0.6; % probability of feature disturb on one sub-graph
param.variationConProb = 0.1; % probability of feature disturb on the whole graph

testPara = [0.1];

for iterTest = 1:nTest
    display([num2str(iterTest)]);
    for iTestPara = 1:length(testPara)
        %%%%%% change this for various test %%%%%
        param.variationConProb = testPara(iTestPara);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        affinity = generateRandomPartitionAffinity(halfNodeCnt,param);
        
        %% try graph cuts
        % raw graph cut
        sumFeature = diag(sum(affinity.featureWeight,1));
        Lap_g = sumFeature - affinity.featureWeight;
        [eigVec, eigD]=eigs(Lap_g,2,'SA');
        if abs(0-eigD(1)) < 1e-6
            Cut_raw = sign(eigVec(:,2));
        else
            Cut_raw = sign(eigVec(:,1));
        end
        
        % balanced graph cut
        if abs(0-eigD(1)) < 1e-6
            Bcut = eigVec(:,2);
        else
            Bcut = eigVec(:,1);
        end
        Median_Bcut = median(Bcut);
        Bcut_raw = zeros(size(Bcut));
        Bcut_raw(find(Bcut>=Median_Bcut))=1;
        Bcut_raw(find(Bcut<Median_Bcut))=-1;
        % [v d] = eigs(Lap_g,1,'LA');
        % tmpLabel = sign(v(:));
        
        
        % [partition] = RayQotCutMatching(Lap,zeros(2*halfNodeCnt,2*halfNodeCnt),1,0,3);
        % partition = real(partition);
        % tmpLabel = sign(partition(:));
        if isDisplay
            tmpLabel = Cut_raw;
%             figure; title('Initial graph cut');
%             subplot(1,2,1); hold on; triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
%             scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),40,'r','MarkerFaceColor','r');
%             scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),40,'b','MarkerFaceColor','b'); hold off;
%             subplot(1,2,2); hold on; triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
%             scatter(affinity.graph(1,find(tmpLabel==1)),affinity.graph(2,find(GClabel==1)),40,'r','MarkerFaceColor','w');
%             scatter(affinity.graph(1,find(tmpLabel==-1)),affinity.graph(2,find(GClabel==-1)),40,'r','MarkerFaceColor','b'); hold off;
        
            figure(6); title('Graph with two partitions'); hold on; axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),40,'r','MarkerFaceColor','r');
            scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),40,'b','MarkerFaceColor','b');
            hold off;
        
            figure(7); title('Initial graph cut result'); hold on; axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,find(tmpLabel==1)),affinity.graph(2,find(tmpLabel==1)),40,'r','MarkerFaceColor','c');
            scatter(affinity.graph(1,find(tmpLabel==-1)),affinity.graph(2,find(tmpLabel==-1)),40,'r','MarkerFaceColor','y');
            hold off;
            
            figure; title('Balanced graph cut result'); hold on; axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,find(Bcut>=Median_Bcut)),affinity.graph(2,find(Bcut>=Median_Bcut)),40,'r','MarkerFaceColor','c');
            scatter(affinity.graph(1,find(Bcut<Median_Bcut)),affinity.graph(2,find(Bcut<Median_Bcut)),40,'r','MarkerFaceColor','y');
            hold off;
        end
        % sumGT = diag(sum(affinity.GT,1));
        % Lap_m = sumGT - affinity.GT;
        % Lap = lambda_1*Lap_g - lambda_2*Lap_m;
        % [v d]=eigs(Lap_g,2,'SA');
        % tmpLabel = sign(v(:,1));
        % y = RayQotCutMatching(affinity.featureWeight, affinity.GT, lambda_1, lambda_2, 1);
        % tmpLabel = sign(y);
        
        % figure(1); title('Initial cut'); hold on;
        % triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
        % scatter(affinity.graph(1,find(tmpLabel==1)),affinity.graph(2,find(tmpLabel==1)),40,'r','MarkerFaceColor','w');
        % scatter(affinity.graph(1,find(tmpLabel==-1)),affinity.graph(2,find(tmpLabel==-1)),40,'r','MarkerFaceColor','b');
        % hold off;
        
        % figure; title('cut using matching cor');
        % subplot(1,2,1); hold on; triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
        % scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),40,'r','MarkerFaceColor','r');
        % scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),40,'b','MarkerFaceColor','b'); hold off;
        % subplot(1,2,2); hold on; triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
        % scatter(affinity.graph(1,find(tmpLabel==1)),affinity.graph(2,find(tmpLabel==1)),40,'r','MarkerFaceColor','w');
        % scatter(affinity.graph(1,find(tmpLabel==-1)),affinity.graph(2,find(tmpLabel==-1)),40,'r','MarkerFaceColor','b'); hold off;
        
        %% CutMatch
        %%%%%%%%%%% 200-50 %%%%%%%%%%%%%%%
        param.lambda_1 = 200; param.lambda_2 = 50;
        [Matching_200_50, Cut_200_50]=CutMatching(affinity.aff, affinity.featureWeight,param);
        Cut_cutmatch_200_50 = sign(Cut_200_50);
        maxWeight = max(Matching_200_50(:));
        [Match_cutmatch_200_50 cost]=hungarian(maxWeight*ones(2*halfNodeCnt,2*halfNodeCnt)-Matching_200_50);
        
        %%%%%%%%%%% 100-20 %%%%%%%%%%%%%%%
        param.lambda_1 = 100; param.lambda_2 = 20;
        [Matching_100_20, Cut_100_20]=CutMatching(affinity.aff, affinity.featureWeight,param);
        Cut_cutmatch_100_20 = sign(Cut_100_20);
        maxWeight = max(Matching_100_20(:));
        [Match_cutmatch_100_20 cost]=hungarian(maxWeight*ones(2*halfNodeCnt,2*halfNodeCnt)-Matching_100_20);
        
        %%%%%%%%%%% 150-30 %%%%%%%%%%%%%%%
        param.lambda_1 = 150; param.lambda_2 = 30;
        [Matching_150_30, Cut_150_30]=CutMatching(affinity.aff, affinity.featureWeight,param);
        Cut_cutmatch_150_30 = sign(Cut_200_50);
        maxWeight = max(Matching_150_30(:));
        [Match_cutmatch_150_30 cost]=hungarian(maxWeight*ones(2*halfNodeCnt,2*halfNodeCnt)-Matching_150_30);
%         
        if isDisplay
            figure; hold on; title('Matching result with CutMatch'); axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),40,'r','MarkerFaceColor','r');
            scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),40,'b','MarkerFaceColor','b');
            for i = 1:2*halfNodeCnt-1
                for j = i:2*halfNodeCnt
                    if Match_cutmatch_200_50(i,j)==1
                        line([affinity.graph(1,i),affinity.graph(1,j)],[affinity.graph(2,i),affinity.graph(2,j)],'LineStyle','--','Color','g','LineWidth',2);
                    end
                end
            end
            hold off;

            figure; title('Cut result with CutMatch'); axis off;
%             subplot(1,2,1); hold on; triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
%             scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),40,'r','MarkerFaceColor','r');
%             scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),40,'b','MarkerFaceColor','b'); hold off;
%             subplot(1,2,2);
            hold on; triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            tmpLabel = Cut_cutmatch_100_20;
            scatter(affinity.graph(1,find(tmpLabel==1)),affinity.graph(2,find(tmpLabel==1)),40,'r','MarkerFaceColor','w');
            scatter(affinity.graph(1,find(tmpLabel==-1)),affinity.graph(2,find(tmpLabel==-1)),40,'r','MarkerFaceColor','k'); hold off;
        end
        %% try graph matching
        X_init = ones([2*halfNodeCnt 2*halfNodeCnt])/(2*halfNodeCnt);
        X_init = BregmanBiStochZero(X_init,50);
        X_init = X_init(:);
        
        X_next = IBGP(X_init, affinity.aff, zeros(2*halfNodeCnt,1), 0, 0.002, 400, 200, 1);
        X_out = reshape(X_next, [2*halfNodeCnt, 2*halfNodeCnt]);
        
        maxWeight = max(X_out(:));
        [Match_ibgp cost] = hungarian(maxWeight*ones(2*halfNodeCnt,2*halfNodeCnt)-X_out);
        
        algpar.iterMax1 = 300;
        X_rrwm_init = 1/(2*halfNodeCnt)*ones(2*halfNodeCnt,2*halfNodeCnt);
        X_rrwm = RRWM(affinity.aff,2*halfNodeCnt, 2*halfNodeCnt,algpar);
        maxWeight2 = max(X_rrwm(:));
        X_rrwm = reshape(X_rrwm,[2*halfNodeCnt,2*halfNodeCnt]);
        [Match_rrwm cost] = hungarian(maxWeight2*ones(2*halfNodeCnt,2*halfNodeCnt)-X_rrwm);
        
        X_ipfp = IPFP(affinity.aff,2*halfNodeCnt,2*halfNodeCnt);
        maxWeightIPFP = max(X_ipfp(:));
        X_ipfp = reshape(X_ipfp,[2*halfNodeCnt,2*halfNodeCnt]);
        [Match_ipfp cost] = hungarian(maxWeightIPFP*ones(2*halfNodeCnt,2*halfNodeCnt)-X_ipfp);
        
        if isDisplay
            figure; hold on; title('Initial graph matching with IBGP'); axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),40,'r','MarkerFaceColor','r');
            scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),40,'b','MarkerFaceColor','b');
            for i = 1:2*halfNodeCnt-1
                for j = i:2*halfNodeCnt
                    if Match_ibgp(i,j)==1
                        line([affinity.graph(1,i),affinity.graph(1,j)],[affinity.graph(2,i),affinity.graph(2,j)],'LineStyle','--','Color','g','LineWidth',2);
                    end
                end
            end
            hold off;

            figure; hold on; title('matching within one graph using RRWM'); axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),25,'r','MarkerFaceColor','r');
            scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),25,'b','MarkerFaceColor','b');
            for i = 1:2*halfNodeCnt-1
                for j = i:2*halfNodeCnt
                    if Match_rrwm(i,j)==1
                        line([affinity.graph(1,i),affinity.graph(1,j)],[affinity.graph(2,i),affinity.graph(2,j)],'LineStyle','--','Color','g','LineWidth',2);
                    end
                end
            end
            hold off;
            
            figure; hold on; title('matching within one graph using IPFP'); axis off;
            triplot(affinity.tri,affinity.graph(1,:), affinity.graph(2,:));
            scatter(affinity.graph(1,1:halfNodeCnt), affinity.graph(2,1:halfNodeCnt),25,'r','MarkerFaceColor','r');
            scatter(affinity.graph(1,halfNodeCnt+1:end), affinity.graph(2,halfNodeCnt+1:end),25,'b','MarkerFaceColor','b');
            for i = 1:2*halfNodeCnt-1
                for j = i:2*halfNodeCnt
                    if Match_ipfp(i,j)==1
                        line([affinity.graph(1,i),affinity.graph(1,j)],[affinity.graph(2,i),affinity.graph(2,j)],'LineStyle','--','Color','g','LineWidth',2);
                    end
                end
            end
            hold off;
        end
        a = 1;
        
        GTmatching = [zeros(halfNodeCnt,halfNodeCnt),eye(halfNodeCnt);eye(halfNodeCnt),zeros(halfNodeCnt,halfNodeCnt)];
        GTcut = [ones(halfNodeCnt,1);-ones(halfNodeCnt,1)];
        diff_matching_cutmatch_200_50 = GTmatching.*Match_cutmatch_200_50;
        diff_matching_cutmatch_150_30 = GTmatching.*Match_cutmatch_150_30;
        diff_matching_cutmatch_100_20 = GTmatching.*Match_cutmatch_100_20;
        diff_matching_rrwm = GTmatching.*Match_rrwm;
        diff_matching_IPFP = GTmatching.*Match_ipfp;
        diff_matching_IBGP = GTmatching.*Match_ibgp;
        diff_cut_num_200_50 = max([sum(Cut_cutmatch_200_50==GTcut),sum(Cut_cutmatch_200_50~=GTcut)]);
        diff_cut_num_150_30 = max([sum(Cut_cutmatch_150_30==GTcut),sum(Cut_cutmatch_150_30~=GTcut)]);
        diff_cut_num_100_20 = max([sum(Cut_cutmatch_100_20==GTcut),sum(Cut_cutmatch_100_20~=GTcut)]);
        diff_cut_numraw = max([sum(Cut_raw==GTcut),sum(Cut_raw~=GTcut)]);
        diff_cut_numBcut = max([sum(Bcut_raw==GTcut),sum(Bcut_raw~=GTcut)]);
        
        acc_matching_200_50{iterTest, iTestPara} = length(find(diff_matching_cutmatch_200_50==1))/(2*halfNodeCnt);
        acc_matching_150_30{iterTest, iTestPara} = length(find(diff_matching_cutmatch_150_30==1))/(2*halfNodeCnt);
        acc_matching_100_20{iterTest, iTestPara} = length(find(diff_matching_cutmatch_100_20==1))/(2*halfNodeCnt);
        acc_cut_200_50{iterTest, iTestPara} = diff_cut_num_200_50/(2*halfNodeCnt);
        acc_cut_150_30{iterTest, iTestPara} = diff_cut_num_150_30/(2*halfNodeCnt);
        acc_cut_100_20{iterTest, iTestPara} = diff_cut_num_100_20/(2*halfNodeCnt);
        acc_matching_rrwm{iterTest, iTestPara} = length(find(diff_matching_rrwm==1))/(2*halfNodeCnt);
        acc_matching_IBGP{iterTest, iTestPara} = length(find(diff_matching_IBGP==1))/(2*halfNodeCnt);
        acc_matching_IPFP{iterTest, iTestPara} = length(find(diff_matching_IPFP==1))/(2*halfNodeCnt);
        acc_cut_raw{iterTest, iTestPara} = diff_cut_numraw/(2*halfNodeCnt);
        acc_Bcut_raw{iterTest, iTestPara} = diff_cut_numBcut/(2*halfNodeCnt);
    end
end

warning on;
for j = 1:length(testPara)
    acc_ave_cutmatch_gm_200_50{j}=0;
    acc_ave_cutmatch_gm_150_30{j}=0;
    acc_ave_cutmatch_gm_100_20{j}=0;
    acc_ave_cutmatch_gc_200_50{j}=0;
    acc_ave_cutmatch_gc_150_30{j}=0;
    acc_ave_cutmatch_gc_100_20{j}=0;
    acc_ave_rrwm_gm{j}=0;
    acc_ave_ipfp_gm{j}=0;
    acc_ave_rawcut_gc{j}=0;
    acc_ave_ibgp_gm{j}=0;
    acc_ave_Bcut_gc{j}=0;
end

for j = 1:length(testPara)
    for i = 1:nTest
        acc_ave_cutmatch_gm_200_50{j} = acc_ave_cutmatch_gm_200_50{j} + acc_matching_200_50{i,j};
        acc_ave_cutmatch_gm_150_30{j} = acc_ave_cutmatch_gm_150_30{j} + acc_matching_150_30{i,j};
        acc_ave_cutmatch_gm_100_20{j} = acc_ave_cutmatch_gm_100_20{j} + acc_matching_100_20{i,j};
        acc_ave_cutmatch_gc_200_50{j} = acc_ave_cutmatch_gc_200_50{j} + acc_cut_200_50{i,j};
        acc_ave_cutmatch_gc_150_30{j} = acc_ave_cutmatch_gc_150_30{j} + acc_cut_150_30{i,j};
        acc_ave_cutmatch_gc_100_20{j} = acc_ave_cutmatch_gc_100_20{j} + acc_cut_100_20{i,j};
        acc_ave_rrwm_gm{j} = acc_ave_rrwm_gm{j} + acc_matching_rrwm{i,j};
        acc_ave_rawcut_gc{j} = acc_ave_rawcut_gc{j} + acc_cut_raw{i,j};
        acc_ave_ibgp_gm{j} = acc_ave_ibgp_gm{j} + acc_matching_IBGP{i,j};
        acc_ave_ipfp_gm{j} = acc_ave_ipfp_gm{j} + acc_matching_IPFP{i,j};
        acc_ave_Bcut_gc{j} = acc_ave_Bcut_gc{j} + acc_Bcut_raw{i,j};
    end
    acc_ave_cutmatch_gm_200_50{j} = acc_ave_cutmatch_gm_200_50{j}/nTest;
    acc_ave_cutmatch_gm_150_30{j} = acc_ave_cutmatch_gm_150_30{j}/nTest;
    acc_ave_cutmatch_gm_100_20{j} = acc_ave_cutmatch_gm_100_20{j}/nTest;
    acc_ave_cutmatch_gc_200_50{j} = acc_ave_cutmatch_gc_200_50{j}/nTest;
    acc_ave_cutmatch_gc_150_30{j} = acc_ave_cutmatch_gc_150_30{j}/nTest;
    acc_ave_cutmatch_gc_100_20{j} = acc_ave_cutmatch_gc_100_20{j}/nTest;
    acc_ave_rrwm_gm{j} = acc_ave_rrwm_gm{j}/nTest;
    acc_ave_rawcut_gc{j} = acc_ave_rawcut_gc{j}/nTest;
    acc_ave_ibgp_gm{j} = acc_ave_ibgp_gm{j}/nTest;
    acc_ave_ipfp_gm{j} = acc_ave_ipfp_gm{j}/nTest;
    acc_ave_Bcut_gc{j} = acc_ave_Bcut_gc{j}/nTest;
end

figure; hold on; title('Matching performance on synthetic data');
p1 = plot(testPara,cell2mat(acc_ave_cutmatch_gm_200_50),'LineWidth',2,'LineStyle','-.','Marker','o','MarkerSize',6,'Color','b');
p2 = plot(testPara,cell2mat(acc_ave_cutmatch_gm_150_30),'LineWidth',2,'LineStyle','-','Marker','s','MarkerSize',6,'Color','y');
p3 = plot(testPara,cell2mat(acc_ave_cutmatch_gm_100_20),'LineWidth',2,'LineStyle','-.','Marker','d','MarkerSize',6,'Color','k');
p4 = plot(testPara,cell2mat(acc_ave_ibgp_gm),'LineWidth',2,'LineStyle','--','Marker','*','MarkerSize',6,'Color','g');
p5 = plot(testPara,cell2mat(acc_ave_rrwm_gm),'LineWidth',2,'LineStyle','-','Marker','p','MarkerSize',6,'Color','r');
p6 = plot(testPara,cell2mat(acc_ave_ipfp_gm),'LineWidth',2,'LineStyle','--','Marker','+','MarkerSize',6,'Color','m');
legend([p1,p2,p3,p4,p5,p6],'CutMatch-200-50','CutMatch-150-30','CutMatch-100-20','IBGP','RRWM','IPFP','Location','Best');
ylabel('\fontname{times new roman}Matching accuracy','FontSize',15);
xlabel('\fontname{times new roman}Graph deformation','FontSize',15);
% text(100,100,'\fontname{times new roman}test \it\sigma = 0.1');
hold off;

figure; hold on; title('Cuts performance on synthetic data');
p1 = plot(testPara,cell2mat(acc_ave_cutmatch_gc_200_50),'LineWidth',2,'LineStyle','-.','Marker','o','MarkerSize',6,'Color','b');
p2 = plot(testPara,cell2mat(acc_ave_cutmatch_gc_150_30),'LineWidth',2,'LineStyle','-','Marker','s','MarkerSize',6,'Color','y');
p3 = plot(testPara,cell2mat(acc_ave_cutmatch_gc_100_20),'LineWidth',2,'LineStyle','-.','Marker','d','MarkerSize',6,'Color','k');
p4 = plot(testPara,cell2mat(acc_ave_rawcut_gc),'LineWidth',2,'LineStyle','--','Marker','*','MarkerSize',6,'Color','g');
p5 = plot(testPara,cell2mat(acc_ave_Bcut_gc),'LineWidth',2,'LineStyle','-.','Marker','v','MarkerSize',6,'Color','c');
legend([p1,p2,p3,p4,p5],'CutMatch-200-50','CutMatch-150-30','CutMatch-100-20','Raw Cuts','Balanced Cuts','Location','Best');
ylabel('\fontname{times new roman}Cuts accuracy','FontSize',15);
xlabel('\fontname{times new roman}Graph deformation','FontSize',15);
hold off;