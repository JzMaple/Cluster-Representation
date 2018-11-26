function plotPairwiseMatching(X,Gt,feat,adjMask,fileName,img)
    nodeCnt = 10;
    figure('Name',fileName,'NumberTitle','off','color','white','Visible','off');
    [height1,~] = size(img{1});
    [height2,~]= size(img{2});
    ratio = height1/height2;
    if ratio>1
        img{1} = imresize(img{1},1/ratio);
        feat{1} = feat{1}/ratio;
    else
        img{2} = imresize(img{2},ratio);
        feat{2} = feat{2}*ratio;
    end
    appImg = appendimages(img{1},img{2});
    feat{2}(:,1) = feat{2}(:,1) + size(img{1},2);
    imshow(appImg); hold on;
    colorCode = makeColorCode(100);
    set(gca,'position',[0 0 1 1]);
    for viewk=1:2
        for j=1:nodeCnt
            plot(feat{viewk}(j,1),feat{viewk}(j,2),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor',colorCode(:,j),'MarkerSize', 10);
            hold on;
        end
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
    xperm = rem(find(X'>0),nodeCnt);
    xperm(xperm==0)=nodeCnt;
    gtperm = rem(find(Gt'>0),nodeCnt);
    gtperm(gtperm==0)=nodeCnt;
    for i=1:nodeCnt
        if xperm(i)~=gtperm(i)
            col = 'r';
        else
            col = 'g';
        end
        plot([ feat{1}(i,1), feat{2}(xperm(i),1) ]...
            ,[ feat{1}(i,2), feat{2}(xperm(i),2) ],...
            '-','LineWidth',2,'MarkerSize',10,...
            'color', col);
    end
    axis tight;
    saveFrame = getframe(gcf);
    saveImg = saveFrame.cdata;
    saveImg = imresize(saveImg,[350,900]);
    imwrite(saveImg,[fileName,'.png']);
    close all;
end

function [ priorColorCode ] = makeColorCode(nCol)

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
end

% im = appendimages(image1, image2)
%
% Return a new image that appends the two images side-by-side.

function im = appendimages(image1, image2)

% Select the image with the fewest rows and fill in enough empty rows
%   to make it the same height as the other image.
rows1 = size(image1,1);
rows2 = size(image2,1);

if (rows1 < rows2)
     image1(rows2,1) = 0;
else
     image2(rows1,1) = 0;
end

% Now append both images side-by-side.
im = [image1 image2];
end