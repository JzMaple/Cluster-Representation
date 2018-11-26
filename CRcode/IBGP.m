function X_next = IBGP(X_init, affinity, y_update, lambda_2, xLearningR, xIterCnt, gIterCnt, mode)
    % IBGP iterative bregman gradeint projection for solving general graph
    % matching
    % solving X_next = argmax_{X}vec(X)'*A*vec(X)+y*y'
    % mode 1: within one graph, 2: standard graph matching
    if ~exist('mode') mode =1; end
    tol_g = 0.003;
    tol_scr = 0.003;
    if lambda_2 == 0
        y_update = zeros(sqrt(size(X_init,1)),1);
    end
    CoupleGradient = -lambda_2*(y_update*y_update');        % calculating n*n matrix gradient on linear part
    nodeCnt = length(y_update);
    X_prev = X_init;
    oneMat = ones(nodeCnt,nodeCnt);
    mask = ones(nodeCnt,nodeCnt)-eye(nodeCnt);
    compRate = ones(nodeCnt,nodeCnt)/xLearningR;
    if mode == 1
        for i = 1:xIterCnt
            GMGradient = (affinity + affinity')*X_prev;  % calculating n^2 vector gradient on quadratic part
            V = reshape(GMGradient,[nodeCnt nodeCnt])+CoupleGradient;
            
            Xmat = reshape(X_prev,[nodeCnt nodeCnt]);
            
            
            V = (V+V')/2;
%            for j = 1:gIterCnt
%                 V_prev = V;
%                 V = V - 1/nodeCnt*V*oneMat - 1/nodeCnt*oneMat*V + 1/(nodeCnt*nodeCnt)*oneMat*V*oneMat;
%                 V = V.*mask;
%                 % V=(V+V')/2;
%                 BoundLower = arrayfun(@max, -Xmat/xLearningR, -Xmat'/xLearningR);
%                 V = arrayfun(@XtruncLower,V,BoundLower);
%                 BoundUpper = arrayfun(@min,compRate-Xmat/xLearningR, compRate-Xmat'/xLearningR);
%                 V = arrayfun(@XtruncUpper,V,BoundUpper);
%                if norm(V - V_prev,'fro') < tol_g
%                    break;
%                end
%             end
            V_prev = V;
            V = mexBistochNormalizeZero(V, 0.003, int32(200), xLearningR, oneMat, Xmat);
            X_next = X_prev + xLearningR*V(:);
            if norm(X_next - X_prev) > tol_scr
                X_prev = X_next;
            else
                break;
            end
        end
    elseif mode == 2
        for i = 1:xIterCnt
            GMGradient = (affinity + affinity')*X_prev;  % calculating n^2 vector gradient on quadratic part
            V = reshape(GMGradient,[nodeCnt nodeCnt])+CoupleGradient;
            
            Xmat = reshape(X_prev,[nodeCnt nodeCnt]);
            % V = (V+V')/2;
            for j = 1:gIterCnt
                V_prev = V;
                V = V - 1/nodeCnt*V*oneMat - 1/nodeCnt*oneMat*V + 1/(nodeCnt*nodeCnt)*oneMat*V*oneMat;
                % V = V.*mask;
                % V=(V+V')/2;
                % BoundLower = arrayfun(@max, -Xmat/xLearningR, -Xmat'/xLearningR);
                V = arrayfun(@XtruncLower,V,-Xmat/xLearningR);
                % BoundUpper = arrayfun(@min,compRate-Xmat/xLearningR, compRate-Xmat'/xLearningR);
                V = arrayfun(@XtruncUpper,V,compRate-Xmat/xLearningR);
                if norm(V - V_prev,'fro') < tol_g
                    break;
                end
            end
            X_next = X_prev + xLearningR*V(:);
            if norm(X_next - X_prev) > tol_scr
                X_prev = X_next;
            else
                break;
            end
        end
    end

end