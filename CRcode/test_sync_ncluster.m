% test the performance of IMGM on synthetic data
clear *; clear -global *; clear all; clc; close all;
global affinity target
varyMinGrhCnt=75; % minimum graph count
varyMaxGrhCnt=75; % maximum graph count
grhTestCnt = 10;  %

setPlotColor;
setObsoleteVariables;
target.config.testType = 'formal'; % to control the test type: formal; massOutlier;

algpar = setPairwiseSolver();
mpmAlgPar = setMPMAlgPar;

target.config.database = 'synthetic';% only synthetic test is allowed here
target.config.Sacle_2D = 0.05;
iterRange = 6;
graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt; %
graphCntList = [graphMinCnt:4:graphMaxCnt]; % sum of all graphs in each test
graphIncreList = [graphMinCnt-22:4:graphMaxCnt-22]; % sum of incremental graphs in each test

if strcmp(target.config.testType,'massOutlier')% outlier test, see Fig.5(b) in the PAMI paper
    % cao cao_c cao_uc cao_pc are not used in massive outlier mode, because
    % no need to enforce consistency by post-step, thus we disable them:
    algNameSepSpace = '                ';
    algSet.algNameSet = {'mpm','rrwm','cao_s','cao_','cao_c_s','cao_c_','cao_uc_s','cao_uc_','cao_pc_s','cao_pc_','mOpt','mSync'};
    algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1,1];
    algSet.algColor = {mpmClr,rrwmClr,caoClr,caoClr,cao_cClr,cao_cClr,cao_ucClr,cao_ucClr,cao_pcClr,cao_pcClr,iccvClr,nipsClr};
    algSet.algLineStyle = {'--','--','--','-','--','-','--','-','--','-','-','-'};
    algSet.algMarker = {'.','.','.','.','.','.','.','.','.','.','.','.'};
    target.config.bGraphMatch = 0;% set to 1 use random graphs, otherwise use random points as set in the MPM code/paper
    target.config.category = 'outlier';% only outlier are supported here
    target.config.inCntType = 'exact';% set 'exact' for "more outlier case", e.g. Fig.5 and Fig.6
    nInlier = 6;target.config.nOutlier = 12;target.config.deform = .05;
    target.config.density = 1;target.config.complete = 1;
    graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
else
    algNameSepSpace = '                    ';
    algSet.algNameSet = {'rrwm','cao_','cao','cao_c_','cao_c','cao_uc_','cao_uc','cao_pc_','cao_pc','mOpt','mSync'};
    algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1];
    algSet.algColor = {rrwmClr,caoClr,caoClr,cao_cClr,cao_cClr,cao_ucClr,cao_ucClr,cao_pcClr,cao_pcClr,iccvClr,nipsClr};
    algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-','-','-'};
    algSet.algMarker = {'.','.','.','.','.','.','.','.','.','.','.'};
    target.config.bGraphMatch = 1;
    target.config.inCntType = 'all';% set 'all' for "only a few outlier case", e.g. Fig.1&2&3&4
    target.config.category = 'outlier';%'deform','outlier','density','complete'
    switch target.config.category
        case 'deform'% same setting with 5th row in Table 1 in the PAMI paper 
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.15;
            target.config.density = .9;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
        case 'outlier'% same setting with 6th row in Table 1 in the PAMI paper 
            nInlier = 6;target.config.nOutlier = 4;target.config.deform = 0.05;
            target.config.density = 1;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
        case 'density'% same setting with 7th row in Table 1 in the PAMI paper 
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.0;
            target.config.density = 0.5;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
        case 'complete'% same setting with 8th row in Table 1 in the PAMI paper 
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.05;
            target.config.density = 1;target.config.complete = 0.1;     
    end
end

graphRange = graphMinCnt:4:graphMaxCnt;
target.config.initConstWeight = .2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constStep = 1.1;% inflate parameter, suggest 1.1-1.2
target.config.constWeightMax = 1;
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3
target.config.edgeAffinityWeight = 1;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:graphMaxCnt;
paraCnt=length(graphRange);
iterCnt = length(iterRange);

[~,rrwmIdx] = ismember('rrwm',algSet.algNameSet);
[~,mpmIdx] = ismember('mpm',algSet.algNameSet);
[~,cao_Idx] = ismember('cao_',algSet.algNameSet);[~,cao_sIdx] = ismember('cao_s',algSet.algNameSet);
[~,caoIdx] = ismember('cao',algSet.algNameSet);
[~,cao_c_Idx] = ismember('cao_c_',algSet.algNameSet);[~,cao_c_sIdx] = ismember('cao_c_s',algSet.algNameSet);
[~,cao_cIdx] = ismember('cao_c',algSet.algNameSet);
[~,cao_uc_Idx] = ismember('cao_uc_',algSet.algNameSet);[~,cao_uc_sIdx] = ismember('cao_uc_s',algSet.algNameSet);
[~,cao_ucIdx] = ismember('cao_uc',algSet.algNameSet);
[~,cao_pc_Idx] = ismember('cao_pc_',algSet.algNameSet);[~,cao_pc_sIdx] = ismember('cao_pc_s',algSet.algNameSet);
[~,cao_pcIdx] = ismember('cao_pc',algSet.algNameSet);
algCnt = length(algSet.algEnable); 
X=cell(algCnt,1);
target.config.nodeCnt = length(target.config.selectNodeMask);
target.config.graphCnt = min(max(graphRange),length(target.config.selectGraphMask{1}));
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;

% paraCnt: iterate over graph #
% algCnt: iterate over algorithms
% testCnt: iterate over tests
timAve = zeros(paraCnt,algCnt,testCnt);timAveFull = zeros(paraCnt,algCnt);
accAve = zeros(paraCnt,algCnt,testCnt);accAveFull = zeros(paraCnt,algCnt);
scrAve = zeros(paraCnt,algCnt,testCnt);scrAveFull = zeros(paraCnt,algCnt);
conPairAve = zeros(paraCnt,algCnt,testCnt);conPairAveFull = zeros(paraCnt,algCnt);
accStd = zeros(paraCnt,algCnt,testCnt);accStdFull = zeros(paraCnt,algCnt);
scrStd = zeros(paraCnt,algCnt,testCnt);scrStdFull = zeros(paraCnt,algCnt);
conPairStd = zeros(paraCnt,algCnt,testCnt);conPairStdFull = zeros(paraCnt,algCnt);

conMatPairGraph = cell(paraCnt,algCnt,testCnt);
accMatPairGraph = cell(paraCnt,algCnt,testCnt);
scrMatPairGraph = cell(paraCnt,algCnt,testCnt);

fidPerf = fopen('results.csv','w');
fprintf(fidPerf, 'testType,bGraphMatch,unaryFeat,edgeFeat,inCntType,database,category,testCnt,iter#,total node#,outlier#,complete,density,deform,graph#,alg#,scale,edgeWeight,initConstWeight,consStep,constWeightMax,iterImmune\n');
fprintf(fidPerf, '%s,%d,%d,%d,%s,%s,%s,%d,%d,%d,%d,%.2f,%.2f,%.2f,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n',...
    target.config.testType,target.config.bGraphMatch,target.config.bUnaryEnable,target.config.bEdgeEnable,target.config.inCntType,target.config.database,target.config.category,...
    testCnt,max(iterRange),nodeCnt,target.config.nOutlier,target.config.complete,target.config.density,target.config.deform,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,target.config.edgeAffinityWeight,...
    target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('testType=%s, bGraphMatch=%d, unaryEnable=%d, edgeEnable=%d, inCntType=%s, database=%s, category=%s, iter#=%d, test#=%d, node#=%d, outlier#=%d,complete=%.2f, density=%.2f, deform=%.2f, graph#=%d, alg#=%d, edgeWeight=%.2f, scale=%.2f, initW=%.2f, stepW=%.2f, maxW=%.2f,iterImmune=%d\n',...
    target.config.testType,target.config.bGraphMatch,target.config.bUnaryEnable,target.config.bEdgeEnable,target.config.inCntType,target.config.database,target.config.category,max(iterRange),testCnt,...
    nodeCnt,target.config.nOutlier,target.config.complete,target.config.density,target.config.deform,graphCnt,sum(algSet.algEnable),...
    target.config.edgeAffinityWeight,target.config.Sacle_2D,target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('\n');fprintf(fidPerf,'\n');

inlierAcc = zeros(paraCnt,testCnt);
% now start to test
estErr = zeros(testCnt,1);

%% performance variables
increGraphCnt = 50;
%%%%%%%%%%%%% 1:raw CAO; 2:incre CAO; 3:IMGM-D; 4:IMGM-R; 5:TIP_Yan;
%%%%%%%%%%%%% 6:mSync;
activeMethod = [1,1,1,1,1,1];
nMethods = 2; % number of all algorithms
accResult = zeros(nMethods,increGraphCnt,grhTestCnt);
conResult = zeros(nMethods,increGraphCnt,grhTestCnt);
scrResult = zeros(nMethods,increGraphCnt,grhTestCnt);
timResult = zeros(nMethods,increGraphCnt,grhTestCnt);
param.n = nInlier + target.config.nOutlier;
param.visualization = 0;
param.iterMax = iterRange;
baseGraphCnt = graphCntList(1) - increGraphCnt;

algNameSepSpace = '                    ';
algSet.algNameSet = {'IMGM-D','IMGM-R','CR-D','CR-R'};
fidPerf = fopen('results.csv','w');
T = clock;
fprintf('time: %d : %d : %f\n', T(4), T(5), T(6));

%% main loop for syn test
for iGraphCnt = 1:length(graphCntList) % for each number of graphs
    baseGraphCnt = graphCntList(iGraphCnt) - increGraphCnt;
    for iTest = 1:grhTestCnt   % for each independent test
        affinity = generateRandomAffinity(nInlier,iTest);
        affinity.GT = repmat(eye(nodeCnt,nodeCnt),graphCnt,graphCnt);%just use identity matrix as ground truth matchings
        
        % rrwm pairwise match, once for all graph pairs
        tStart = tic;% compute all pairwise matchings at one time -> rawMat
        rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,iTest);% generate matchings by pairwise matching solver
        rawTotalTime = toc(tStart);
        
        switch target.config.inCntType
            case 'exact' % already known, used in Fig.5 and top two rows in Fig.6
                target.config.inCnt = nodeCnt - target.config.nOutlier;
            case 'all' % in case of few outliers, used in Fig.1,2,3,4
                target.config.inCnt = nodeCnt;
            case 'spec' % specified by user, used in the bottom row of Fig.6
                target.config.inCnt = specNodeCnt;
        end
        scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
        scrDenomMatInCntGT = cal_pair_graph_inlier_score(affinity.GT,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    
        
        IMGMcount = increGraphCnt;
        scrDenomCurrent = max(max(scrDenomMatInCnt(1:end,1:end)));
        
        for i = 1:IMGMcount
            if i == 1
                baseMat = CAO(rawMat(1:end-nodeCnt*IMGMcount,1:end-nodeCnt*IMGMcount),nodeCnt, baseGraphCnt, iterRange,scrDenomCurrent, 'exact',1);
            end
                     
            %%%%%%%%%%%%%%%%%%%%% CRnew %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if i == 1
                CRnewPrevMatching = baseMat; % DPP 
                CRnewPrevMatchingRnd = baseMat; % random
            end
            
            CRnewMatTmp = rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i));
            CRnewMatTmp(1:end-nodeCnt,1:end-nodeCnt)=CRnewPrevMatching;
            
            CRnewMatRndTmp = rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i));
            CRnewMatRndTmp(1:end-nodeCnt,1:end-nodeCnt)=CRnewPrevMatchingRnd;
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(CRnewMatTmp,affinity.GT(1:nodeCnt*(baseGraphCnt+i),1:nodeCnt*(baseGraphCnt+i)),nodeCnt,baseGraphCnt+i,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(CRnewMatTmp,nodeCnt,baseGraphCnt+i,0);
            
            scrDenomMatInCntRndTmp = cal_pair_graph_inlier_score(CRnewMatRndTmp,affinity.GT(1:nodeCnt*(baseGraphCnt+i),1:nodeCnt*(baseGraphCnt+i)),nodeCnt,baseGraphCnt+i,nodeCnt);
            conDenomMatInCntRndTmp = cal_pair_graph_consistency(CRnewMatRndTmp,nodeCnt,baseGraphCnt+i,0);
            
            %%%%%%%%%%%%%%%%%% AP inital arguments %%%%%%%%%%%%%%%%%%
            sigma = 0;
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            simRnd = (1-sigma)*scrDenomMatInCntRndTmp + sigma*conDenomMatInCntRndTmp;
            param.N = baseGraphCnt + i - 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%% test DPP %%%%%%%%%%%%%%%%%
            param.method = 1;
            tStart = tic;
            CRnewIncreMatching = CRnew(simAP, CRnewMatTmp, param);
            tEnd = toc(tStart);
            CRnewPrevMatching = CRnewIncreMatching;
            
            accImgmU = cal_pair_graph_accuracy(CRnewIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
            scrImgmU = cal_pair_graph_score(CRnewIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
            conImgmU = cal_pair_graph_consistency(CRnewIncreMatching,nodeCnt,baseGraphCnt+i,0);
            
            accResult(1,i,iTest)=mean(accImgmU(:));
            scrResult(1,i,iTest)=mean(scrImgmU(:));
            conResult(1,i,iTest)=mean(conImgmU(:));
            timResult(1,i,iTest)=tEnd;
            
            %%%%%%%%%%%%%% test Random %%%%%%%%%%%%%%
            param.method = 3;
            tStart = tic;
            CRnewIncreMatchingRnd = CRnew(simRnd, CRnewMatRndTmp, param);
            tEnd = toc(tStart);
            CRnewPrevMatchingRnd = CRnewIncreMatchingRnd;
            
            accImgmR = cal_pair_graph_accuracy(CRnewIncreMatchingRnd,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
            scrImgmR = cal_pair_graph_score(CRnewIncreMatchingRnd,affinity.GT,nodeCnt,baseGraphCnt+i);
            conImgmR = cal_pair_graph_consistency(CRnewIncreMatchingRnd,nodeCnt,baseGraphCnt+i,0);
            
            accResult(2,i,iTest)=mean(accImgmR(:));
            scrResult(2,i,iTest)=mean(scrImgmR(:));
            conResult(2,i,iTest)=mean(conImgmR(:));
            timResult(2,i,iTest)=tEnd;
        end
        fprintf('--------------------------------------------------------------test %02d performance-------------------------------------------------------------------\n',iTest);
        fprintf(fidPerf,'test%02d\n',iTest);
        fprintf(algNameSepSpace); fprintf(fidPerf,',,');
        for algk=1:nMethods
            fprintf([algSet.algNameSet{algk},algNameSepSpace]);
            fprintf(fidPerf,[algSet.algNameSet{algk},',,,,']);
        end
        fprintf('\n');fprintf(fidPerf,'\n');
        fprintf('grh# itr#  ');fprintf(fidPerf,'grh#, itr#');
        for i = 1:nMethods
            fprintf(' acc   scr   con   tim   ');
            fprintf(fidPerf,', acc,  score, consis, time');
        end
        fprintf('\n');fprintf(fidPerf,'\n');
        for i = 5:5:50
            fprintf(' %02d,  %02d ',i,iTest);fprintf(fidPerf,' %02d,  %02d ',i,iTest);
            for algk=1:nMethods
                fprintf('| %.3f %.3f %.3f %.3f',accResult(algk,i,iTest),scrResult(algk,i,iTest),conResult(algk,i,iTest),timResult(algk,i,iTest));
                fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f',accResult(algk,i,iTest),scrResult(algk,i,iTest),conResult(algk,i,iTest),timResult(algk,i,iTest));
            end
            fprintf('\n');fprintf(fidPerf,'\n');
        end
        T = clock;
        fprintf('time: %d : %d : %f\n', T(4), T(5), T(6));
    end
end

fprintf('--------------------------------------------------------------overall performance-------------------------------------------------------------------\n');
fprintf(fidPerf,'overall mean\n');
fprintf(algNameSepSpace); fprintf(fidPerf,',,');
for algk=1:nMethods
    fprintf([algSet.algNameSet{algk},algNameSepSpace]);
    fprintf(fidPerf,[algSet.algNameSet{algk},',,,,']);
end
fprintf('\n');fprintf(fidPerf,'\n');
fprintf('grh# itr#  ');fprintf(fidPerf,'grh#, itr#');
for i = 1:nMethods
    fprintf(' acc   scr   con   tim   ');
    fprintf(fidPerf,', acc,  score, consis, time');
end
fprintf('\n');fprintf(fidPerf,'\n');
for i = 1:50
    fprintf(' %02d,  all ',i);fprintf(fidPerf,' %02d,  all',i);
    for algk=1:nMethods
        acc = mean(accResult(algk,i,:));
        scr = mean(scrResult(algk,i,:));        
        con = mean(conResult(algk,i,:));
        tim = mean(timResult(algk,i,:));
        fprintf('| %.3f %.3f %.3f %.3f',acc,scr,con,tim);
        fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f',acc,scr,con,tim);
    end
    fprintf('\n');fprintf(fidPerf,'\n');
    if mod(i + 25, 15) == 0 
        fprintf('--------------------------------------------------------------\n')
    end
end