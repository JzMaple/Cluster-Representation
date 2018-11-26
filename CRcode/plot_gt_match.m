function plot_gt_match(target,titleStr)
imgInput = appendimages(target{1}.img, target{2}.img );
height3 = size(target{3}.img,1);
width3 = size(target{3}.img,2);
upperheight = size(imgInput,1);
lowwidth = max([size(imgInput,2),width3]);
img3 = zeros(height3,lowwidth,size(imgInput,3));
start = max([0,floor((lowwidth-width3)/2)]);
img3(1:height3,start+1:start+width3,:)=target{3}.img;
imgInput = appendimages(imgInput,img3,'v');
imgInput = double(imgInput)./255;
acc{1} = 1;acc{2} = 1;acc{3} = 1;
acc_total = 1;
hd = figure('Name',[titleStr,...
    ' ', num2str(acc_total),' acc12=', num2str(acc{1}),' acc13=', num2str(acc{2}),' acc23=', num2str(acc{3})]...
    ,'NumberTitle','off','color','white');
hold off;iptsetpref('ImshowBorder','tight');clear clf;
imshow(imgInput); hold on;
feat{1} = target{1}.point;
feat{2} = target{2}.point;
feat{3} = target{3}.point;
feat{2}(:,1) = feat{2}(:,1) + size(target{1}.img,2);
feat{3}(:,1) = feat{3}(:,1) + start;
feat{3}(:,2) = feat{3}(:,2) + upperheight;

colorCode = makeColorCode(100);
for viewk=1:3
    for j=1:length(feat{viewk})
        plot(feat{viewk}(j,1),feat{viewk}(j,2),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor',colorCode(:,j),'MarkerSize', 10);
            text(feat{viewk}(j,1),feat{viewk}(j,2),num2str(j));hold on;
    end
end

for pairk=1:3
    if strcmp(titleStr,'ground truth')
        if pairk==1,x=1;y=2;X = target{1}.gtPerm*target{2}.gtPerm';GT = target{1}.gtPerm*target{2}.gtPerm';end
        if pairk==2,x=1;y=3;X = target{1}.gtPerm*target{3}.gtPerm';GT = target{1}.gtPerm*target{3}.gtPerm';end
        if pairk==3,x=2;y=3;X = target{2}.gtPerm*target{3}.gtPerm';GT = target{2}.gtPerm*target{3}.gtPerm';end
    else
        if pairk==1,x=1;y=2;X = target{1}.X;GT = target{1}.gtPerm*target{2}.gtPerm';end
        if pairk==2,x=1;y=3;X = target{2}.X;GT = target{1}.gtPerm*target{3}.gtPerm';end
        if pairk==3,x=2;y=3;X = target{3}.X;GT = target{2}.gtPerm*target{3}.gtPerm';end
    end
%     GT = eye(length(X));
    
    width = size(X,2);
    xidx2 = find(X'>0);
    xidx = rem(xidx2,width);
    xidx(xidx==0)=width;
    gtidx2 = find(GT'>0);
    gtidx = rem(gtidx2,width);
    gtidx(gtidx==0)=width;
    xlen = length(xidx);
    gtlen = length(gtidx);
    assert(gtlen==xlen);
    
    for i = 1:xlen 
        if xidx(i)~= gtidx(i)%mismatch in red
            col = 'r';
        else %right match in yellow
            col = 'y';
        end 
        plot([ feat{x}(i,1), feat{y}(xidx(i),1) ]...
            ,[ feat{x}(i,2), feat{y}(xidx(i),2) ],...
                    '-','LineWidth',2,'MarkerSize',10,...
                    'color', col);
    end
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

nMore = nCol - size(priorColorCode,1);
if nMore > 0 
    priorColorCode(size(priorColorCode,1)+1:nCol,:) = rand(nMore, 3);
end

priorColorCode = priorColorCode';

end