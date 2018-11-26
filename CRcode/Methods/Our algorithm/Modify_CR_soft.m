function [ X ] = CR_soft( M, group1, group2, M1, varargin )
%Graph Matching by CR combine soft assign methods
%  by Yu Tian, 2011

%% Default Parameters
param = struct( ...
    'c', 0.2, ...                   % prob. for walk or reweighted jump?
    'amp_max', 30, ...              % maximum value for amplification procedure
    'iterMax', 300, ...             % maximum value for amplification procedure
    'thresConvergence', 1e-25, ...  % convergence threshold for random walks
    'tolC', 1e-3 ...                % convergence threshold for the Sinkhorn method
);
param = parseargs(param, varargin{:});

%% parameter structure -> parameter value
strField = fieldnames(param);
for i = 1:length(strField), eval([strField{i} '=param.' strField{i} ';']); end

% get groups for bistochastic normalization
[idx1 ID1] = make_groups(group1);
[idx2 ID2] = make_groups(group2);

if ID1(end) < ID2(end)
    [idx1 ID1 idx2 ID2 dumVal dumSize] = make_groups_slack(idx1, ID1, idx2, ID2);
    dumDim = 1;
elseif ID1(end) > ID2(end)
    [idx2 ID2 idx1 ID1 dumVal dumSize] = make_groups_slack(idx2, ID2, idx1, ID1);
    dumDim = 2;
else
    dumDim = 0; dumVal = 0; dumSize = 0;
end
idx1 = idx1-1; idx2 = idx2-1;

%% eliminate conflicting elements to prevent conflicting walks
% conf1 = zeros(size(M));
% conf2 = zeros(size(M));
% for i = 1:size(group1,2)
%     idx = find(group1(:,i));
%     for j = 1:length(idx)
%         for k = 1:length(idx)
%             conf1(idx(j),idx(k)) = 1;
%         end
%     end
% end
% for i = 1:size(group2,2)
%     idx = find(group2(:,i));
%     for j = 1:length(idx)
%         for k = 1:length(idx)
%             conf1(idx(j),idx(k)) = 1;
%         end
%     end
% end
% conf = conf1 | conf2;
% M = M.*(~conf);
%%
nMatch = length(M);
d = max(max(M));
lambda = d/3000*sqrt(nMatch);
if 0
    I = eye(nMatch);
else
    I = M1;
end
%M = M-2*I;

prev_assign = ones(nMatch,1)/nMatch; 
%prev_assign = reshape(prev_assign, nSize, 1);
prev_assign2 = prev_assign;
cur_assign = prev_assign;

final_assign = cur_assign; % 这两句后加的
final_assign2 = cur_assign;
cur_score = 0;

bCont = 0;
iter = 0;
amp_max = 30;
iterMax = 500;
thresConvergence = 1e-25;
tolC = 1e-3;


%%
% 要实现的方法：  1.XtMX-->(匈牙利)YtMZ 然后考虑 1.(M+Lambda*I)Y-->...直到收敛
%                                              2.(M+Lambda*I)Z-->...直到收敛 
%                2.XtMX-->1.（匈牙利一步）YtMX 然后考虑 (M+Lambda*I)Y-->...直到收敛
%                         2.（Soft Assign）YtMX 然后考虑 (M+Lambda*I)Y-->...直到收敛
%                3.Lambda 可以 1.固定步长
%                              2. 自适应+固定步长
%                4.可以考虑Xt(M+Lambda*D)X，D衡量点之间的近似度

%% start main iteration
while 0 == bCont && iter < iterMax
    iter = iter+1;
    cur_result = M*cur_assign;
    sumCurResult = sum(cur_result); % normalization of sum 1
    if sumCurResult>0, cur_result = cur_result./sumCurResult; end
    %cur_assign_mat = postHungarian( problem, cur_result);
    amp_value = amp_max/max(cur_result);
    %amp_value = 10*iter;
    cur_assign = exp( amp_value*cur_result);
    
    X_slack = [cur_assign; dumVal*ones(dumSize,1)];
    X_slack = mexBistocNormalize_match_slack(X_slack, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC, dumDim, dumVal, int32(1000));
    cur_assign = X_slack(1:nMatch);
    
    sumCurAssign = sum(cur_assign);
    if sumCurAssign>0, cur_assign = cur_assign./sumCurAssign; end
    if cur_assign'*M*cur_assign >= cur_score
        final_assign = cur_assign;
        final_assign2 = prev_assign;
    end
    
    if sum((cur_assign-prev_assign).^2) < thresConvergence % 说明得到了不变的迭代解
        bCont = 1;
    elseif sum((cur_assign-prev_assign2).^2) < thresConvergence % 说明进入了 cur-prev循环
        prev_assign = final_assign;
        prev_assign2 = final_assign2;
        M = M + lambda*I; % 通过增加对角 试图跳出循环
    else
        prev_assign2 = prev_assign;
        prev_assign = cur_assign;
    end
end
%X = reshape(final_assign, n1, n2);% cur_assign_mat; % 如果是个好的解 早就存到cur_assign_mat了 返回的是矩阵形式 不是列向量形式
X = cur_assign;

end

function [X,Xslack]=bistocNormalize_slack(X,tolC)
[n1,n2]=size(X);
if n1~=n2
    Xslack=X;
    if n1>n2
        Xslack(:,n2+1:n1)=1;
    else
        Xslack(n1+1:n2,:)=1;
    end
    Xslack = bistocNormalize(Xslack,tolC,1000);
    X=Xslack(1:n1,1:n2);
else
    Xslack=X;
    X = bistocNormalize(X,tolC,1000);
end
end
