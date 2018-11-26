% Bregman projection onto a matrix X under constraints:
% X is bi-stochastic, X_ii = 0, X = X^T
function X = BregmanBiStochCut(X_raw, y, iterMax)
    y = sign(y);
    X_raw = (X_raw+X_raw')/2;
    mask = ones(size(X_raw));
    for i = 1:length(y)
        for j = 1:length(y)
            if y(i) == y(j)
                mask(i,j)=0;
            end
        end
    end
    n = size(X_raw,1);
    % oneVec = ones(n,1);
    oneMat = ones(n,n);
    X = X_raw.*mask;
    mask = ones(length(y),length(y))-eye(length(y));
    for i = 1:iterMax
        % X = X + 1/n*(oneMat-X*oneMat-oneMat*X')+1/(2*n*n)*(oneMat*X*oneMat+oneMat*X'*oneMat);
        X = X + 1/n*(oneMat-X*oneMat-oneMat*X)+1/(n*n)*(oneMat*X*oneMat);
        X(X<0) = 0;
        X = X.*mask;
    end
end