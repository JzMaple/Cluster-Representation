% Bregman projection onto a matrix X under constraints:
% X is bi-stochastic, X_ii = 0, X = X^T
function X = BregmanBiStochZero(X_raw, iterMax)
    mask = ones(size(X_raw)) - eye(size(X_raw,1));
    n = size(X_raw,1);
    oneVec = ones(n,1);
    oneMat = ones(n,n);
    X = X_raw;
    for i = 1:iterMax
        X = X + 1/n*(oneMat-X*oneMat-oneMat*X')+1/(2*n*n)*(oneMat*X*oneMat+oneMat*X'*oneMat);
        X(X<0) = 0;
        X = X.*mask;
    end
end