function [ X ] = GAGM( M, group1, group2, varargin )
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

%%
nMatch = length(M);

prev_assign = ones(nMatch,1)/nMatch; 
%prev_assign = reshape(prev_assign, nSize, 1);
prev_assign2 = prev_assign;
cur_assign = prev_assign;

final_assign = cur_assign; % �������ӵ�

bCont = 0;
iter = 0;
iterMax = 500;
thresConvergence = 1e-25;
tolC = 1e-3;

%% start main iteration
while 0 == bCont && iter < iterMax
    iter = iter+1;
    cur_result = M*cur_assign;
    sumCurResult = sum(cur_result); % normalization of sum 1
    if sumCurResult>0, cur_result = cur_result./sumCurResult; end
    %amp_value = 1000;
    amp_value = 20/max(cur_result);
    cur_assign = exp( amp_value*cur_result);
    
    X_slack = [cur_assign; dumVal*ones(dumSize,1)];
    X_slack = mexBistocNormalize_match_slack(X_slack, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC, dumDim, dumVal, int32(1000));
    cur_assign = X_slack(1:nMatch);
    
    sumCurAssign = sum(cur_assign);
    if sumCurAssign>0, cur_assign = cur_assign./sumCurAssign; end

    if sum((cur_assign-prev_assign).^2) < thresConvergence % ˵���õ��˲���ĵ�����
        final_assign = cur_assign;
        bCont = 1;

    elseif sum((cur_assign-prev_assign2).^2) < thresConvergence % ˵�������� cur-prevѭ��
        if prev_assign'*M*prev_assign > cur_assign'*M*cur_assign % YtMY > XtMX ����ط���M�Ϳ��� ���ó�problem.M ��Ϊ�Ͳ��������
            final_assign = prev_assign;
        else % XtMX > YtMY
            final_assign = cur_assign;
        end
        bCont = 1;
    else
        prev_assign2 = prev_assign;
        prev_assign = cur_assign;
    end
end

%X = reshape(final_assign, n1, n2);% cur_assign_mat; % ����Ǹ��õĽ� ��ʹ浽cur_assign_mat�� ���ص��Ǿ�����ʽ ������������ʽ
X = final_assign;

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
