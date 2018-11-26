a1 = 5; a2 = 15;
W = affinity.K{a1,a2};
W1 = sum(W,1); W2 = sum(W,2);
[nRow nCol] = size(W);
nRow = sqrt(nRow); nCol = sqrt(nCol);
algpar.iterMax1 = 100;

A = zeros(2*nRow, nRow*nRow);
for i = 1:nRow
    for j = 1:nCol
        ind_tmp = sub2ind([nRow nCol],i,j);
        A(i,ind_tmp) = 1;
    end
end

for j = 1:nCol
    for i = 1:nRow
        ind_tmp = sub2ind([nRow nCol],j,i);
        A(i+nRow,ind_tmp) = 1;
    end
end

W_sum = -(W1+W2');
b = ones(1,2*nRow);
lb = zeros(1,nRow*nRow);
ub = ones(1,nRow*nRow);
tic;
x = linprog(W_sum,A,b,[],[],lb,ub);
toc
x_slo = reshape(x,[nRow nCol]);

tic;
x_rr = RRWM(W,affinity.nP{a1},affinity.nP{a2},algpar);
max_weight = max(x_rr(:));
x_rr = reshape(x_rr,[nRow nCol]);
x_rr = max_weight*ones(nRow,nCol) - x_rr;
x_rrwm = hungarian(x_rr);
toc

figure; hold on;
subplot(1,2,1); imshow(x_slo,[]);
subplot(1,2,2); imshow(x_rrwm,[]);
a = 1;