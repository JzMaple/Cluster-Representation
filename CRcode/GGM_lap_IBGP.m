function x = GGM_lap(affinity, param)
    param.delta = 0.05;
    param.iterMax = 50;
    param.iterStep = 0.01;
    param.iterMaxBi = 50;
    [nSize,~] = size(affinity);
    nNode = sqrt(nSize);
    
    x_init = 1/nNode*ones(nNode,nNode);
    x_prev = x_init(:);
    
    oneMat = ones(nNode,nNode);
    
    for i=1:iterMax
        %%% calculating gradient of approximation func %%%
        lap_temp = (exp((x_prev-1)/param.delta)).^2;
        x_grad = 1/param.delta*(affinity*lap_temp);
        X_mat = reshape(x_grad,[nNode,nNode]);
        V = X_mat;
        % x_cur = x_prev + param.iterStep*x_grad;
        % v_prev = reshape(x_cur,[nNode,nNode]);
        %%% Bistoch %%%
        for j = 1:param.iterMaxBi
            v_cur = v_prev + 1/nNode*oneMat-1/nNode*v_prev*oneMat - 1/nNode*oneMat*v_prev + 1/(nNode*nNode)*oneMat*v_prev*oneMat;
            v_cur(v_cur<0) = 0;
            v_prev = v_cur;
        end
        x_cur = v_cur(:);
        if i == param.iterMax
            break;
        end
        x_prev = x_cur;
    end
    
    x = x_cur;
end