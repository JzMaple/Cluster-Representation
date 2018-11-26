function plotPairMatching(X,GT,algName)
global target affinity
methodCnt = length(X);
% bDisplayPause = target.display.bDisplayPause;
bDisplayAdj = target.display.bDisplayAdj;
% bSaveimg = target.display.bSaveimg;
outputPath = target.display.outputPath;
database = target.config.database;
category = target.config.category;
metaAlg = target.config.metaAlg;
adj = affinity.adj;
graphCnt = target.config.graphCnt;
nOutlier = target.config.nOutlier;
nodeCnt = target.config.nodeCnt;
% bPermute = target.config.bPermute;

for i=1:methodCnt
    switch algName{i}{1}(1:end-2)
        case 'r'%raw 被截断了
            rawIdx = i;
        case 'IC'%ICML 被截断了
            icmlIdx = i;
        case 'alg1ms'
            alg1msIdx = i;
        case 'IC'%iccv 被截断了
            iccvIdx = i;
         case 'NI'%iccv 被截断了
            nipsIdx = i;
        case 'alg2m'
            alg2mIdx = i;
        case 'alg3mg'
            alg3mgIdx = i;
        case 'alg3mp'
            alg3mpIdx = i;
        case 'mpm'
            mpmIdx = i;
    end
end

if ~exist(outputPath,'dir'),mkdir(outputPath);end

filename = cell(graphCnt,1);
for i=1:graphCnt
    filename{i} = target.config.filename{i};
end
close;
colorCode = makeColorCode(100);
img = cell(2,1);%pairwise display
feat= cell(2,1);
adjMask = cell(2,1);
P = cell(methodCnt);
acc = zeros(methodCnt,1);
height = zeros(2,1);
width = zeros(2,1);
for x=1:graphCnt
    for y=x+1:graphCnt
        
        posx = strmatch(filename{x},target.display.plotGraphList);
        posy = strmatch(filename{y},target.display.plotGraphList);
        if length(posx)>2||length(posy)>2, continue;end
        xscope = (x-1)*nodeCnt+1:x*nodeCnt;
        yscope = (y-1)*nodeCnt+1:y*nodeCnt;
        Gt = GT(xscope,yscope);
        for mk=1:methodCnt
            P{mk} = X{mk}(xscope,yscope);
            acc(mk) = cal_acc(P{mk},nOutlier,Gt);
        end
        [sortAcc idx] = sort(acc);%从小到大
        
%         if sortAcc(end) - sortAcc(end-1)<0.3||max(acc)~=acc(tipMk)||acc(tipMk)<0.5||min(acc)>0.8,continue;end%如果不是最佳精度，就不保存
        
        
        switch category
            case 'duck'%RRWM ok
                if acc(alg1msIdx)<acc(iccvIdx)||acc(alg1msIdx)<acc(mpmIdx),continue;end
                if max(acc)~=acc(alg2mIdx),continue;end%如果不是最佳精度，就不保存
            case 'motorbike'%GAGM ok
                if acc(alg1msIdx)<acc(iccvIdx)||acc(alg1msIdx)<acc(mpmIdx),continue;end
                if max(acc)~=acc(alg2mIdx),continue;end%如果不是最佳精度，就不保存
            case 'car'%FGM ok
                if acc(alg1msIdx)<acc(iccvIdx)||acc(alg1msIdx)<acc(mpmIdx),continue;end
                if max(acc)~=acc(alg2mIdx),continue;end%如果不是最佳精度，就不保存
            case 'hotel'%FGM ok
                if acc(alg1msIdx)<acc(nipsIdx)||acc(alg1msIdx)>=acc(alg2mIdx),continue;end
                if max(acc)~=acc(alg2mIdx),continue;end%如果不是最佳精度，就不保存
            case 'house'%RRWM ok
                if acc(alg1msIdx)<acc(iccvIdx)||acc(alg1msIdx)>=acc(alg2mIdx),continue;end
            case {'volvoc70','houseblack'}%IPFP ok
                if acc(alg1msIdx)<acc(icmlIdx)||acc(alg1msIdx)>=acc(alg2mIdx),continue;end
        end
        
        target.display.plotGraphList = [target.display.plotGraphList,filename{x},filename{y}];
        
        img{1} = target.data{x}.img;
        img{2} = target.data{y}.img;
           
        [height(x) width(x)] = size(target.data{x}.img);
        [height(y) width(y)]= size(target.data{y}.img);
        ratio = height(x)/height(y);
        feat{1} = target.data{x}.point;
        feat{2} = target.data{y}.point;
        adjMask{1} = adj{x};
        adjMask{2} = adj{y};
        if ratio>1
            img{1} = imresize(img{1},1/ratio);
            feat{1} = feat{1}/ratio;
        else
            img{2} = imresize(img{2},ratio);
            feat{2} = feat{2}*ratio;
        end
        appImg = appendimages(img{1},img{2});
        feat{2}(:,1) = feat{2}(:,1) + size(img{1},2);%第一个图的宽

        for k=1:methodCnt
            fileName = [outputPath,'/',database,'_',category,'_',filename{x},'_',filename{y},'_',algName{k}{1},'_',metaAlg,'_acc=',num2str(acc(k))];
            if exist(fileName,'file'),break;end
            
            
        figure('Name',fileName,'NumberTitle','off','color','white','Visible','off');
        
        imshow(appImg); hold on;
        set(gca,'position',[0 0 1 1]);
        for viewk=1:2%画待匹配两个图的点
            for j=1:nodeCnt-nOutlier%先画内点
                plot(feat{viewk}(j,1),feat{viewk}(j,2),'o','MarkerEdgeColor','k',...
                        'MarkerFaceColor',colorCode(:,j),'MarkerSize', 10);
%                         text(feat{viewk}(j,1),feat{viewk}(j,2),num2str(j));
                hold on;
            end
            for j=nodeCnt-nOutlier+1:nodeCnt%再画外点
                plot(feat{viewk}(j,1),feat{viewk}(j,2),'o','MarkerEdgeColor','k',...
                        'MarkerFaceColor',[1 1 1],'MarkerSize', 10);
%                         text(feat{viewk}(j,1),feat{viewk}(j,2),num2str(j));
                hold on;
            end
            %再画adj
            if bDisplayAdj
                for i=1:nodeCnt
                    for j=i+1:nodeCnt
                        if adjMask{viewk}(i,j)
                            plot([ feat{viewk}(i,1), feat{viewk}(j,1)]...
                            ,[ feat{viewk}(i,2), feat{viewk}(j,2)],...
                            '-','LineWidth',2,'MarkerSize',10,...
                            'color', 'y');
                        end
                    end
                end
                hold on;
            end
        end

        %下面画匹配结果
            xperm{x,y,k} = rem(find(P{k}'>0),nodeCnt);
            xperm{x,y,k}(xperm{x,y,k}==0)=nodeCnt;
            gtperm{x,y,k} = rem(find(Gt'>0),nodeCnt);
            gtperm{x,y,k}(gtperm{x,y,k}==0)=nodeCnt;

            for i=1:nodeCnt-nOutlier
                if xperm{x,y,k}(i)~=gtperm{x,y,k}(i)
                     col = 'r';
                else %正确匹配
                    col = 'g';
                end 
                plot([ feat{1}(i,1), feat{2}(xperm{x,y,k}(i),1) ]...
                    ,[ feat{1}(i,2), feat{2}(xperm{x,y,k}(i),2) ],...
                            '-','LineWidth',2,'MarkerSize',10,...
                            'color', col);
            end
%             print(gcf,'-dpng',[fileName,'.png']) 
            axis tight;
            saveFrame = getframe(gcf);
            saveImg = saveFrame.cdata;
            saveImg = imresize(saveImg,[350,900]);
            imwrite(saveImg,[fileName,'.png']);
            disp([fileName,'.png']);
        end%for methodk
         close all;
    end
end


function [ priorColorCode ] = makeColorCode( nCol )

priorColorCode(1,:) = [ 1 0 0 ]; 
priorColorCode(2,:) = [ 0 1 0 ]; 
priorColorCode(3,:) = [ 0 0 1 ]; 
priorColorCode(4,:) = [ 0 1 1 ]; 
priorColorCode(5,:) = [ 1 0 1 ]; 
priorColorCode(6,:) = [ 1 1 0 ]; 
priorColorCode(7,:) = [ 1 0.5 0 ]; 
priorColorCode(8,:) = [ 1 0 0.5 ]; 
priorColorCode(9,:) = [ 1 0.5 0.5 ]; 
priorColorCode(10,:) = [ 0.5 1 0 ]; 
priorColorCode(11,:) = [ 0 1 0.5 ]; 
priorColorCode(12,:) = [ 0.5 1 0.5 ]; 
priorColorCode(13,:) = [ 0.1 .1 0.1 ]; 
priorColorCode(14,:) = [ 0.9 .9 0.9 ]; 
priorColorCode(15,:) = [ 0.5 0 1 ];
priorColorCode(17,:) = [ 0.5 0.5 1 ];
priorColorCode(18,:) = [ 1 0 0.5 ];
priorColorCode(19,:) = [ 1 0.5 0 ];
priorColorCode(20,:) = [ 0.1 0.5 1 ];
% nMore = nCol - size(priorColorCode,1);
% if nMore > 0 
%     priorColorCode(size(priorColorCode,1)+1:nCol,:) = priorColorCode(nMore, 3);
% end
while length(priorColorCode)<nCol
    priorColorCode = [priorColorCode;priorColorCode];
end
priorColorCode = priorColorCode';