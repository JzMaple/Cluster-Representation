clear all; clc; close all;
testMode = 2;
nodeCnt = 10;
m = [0.1 0.2 0.3 0.4];

aff = rand(100,100);

%% mode 1 and 2: concave function; mode 3: arbitrary function;
if testMode==1
    affRnd = -aff'*aff; % negative semidefinit and symmetric
elseif testMode == 2
    d = diag(rand([100 1]));
    affRnd = -aff*d/aff; % negative semidefinite but asymmetric
elseif testMode == 3
    affRnd = aff; % arbitrary function
end

%% from seed point 1
X_init = generateRandBStoch(10,m); % generating seed 1
X_init_1 = X_init;
X_init = (X_init+X_init')/2;

X_init = X_init(:);
X_init_1=X_init_1(:);

X_next = IBGP(X_init, affRnd, zeros(10,1), 0, 0.002, 400, 200, 2);
X_mat1 = reshape(X_next,[10 10]);
figure(1);
subplot(1,2,1), imshow(X_mat1,[]), title('result from seed 1');
subplot(1,2,2), imshow(reshape(X_init,[10 10]),[]), title('seed 1');
display(['affinity score IBGP seed 1: ',num2str(X_next'*affRnd*X_next)]);

%% from seed point 2
X_init = generateRandBStoch(10,m); % generating seed 2
X_init_2 = X_init;
X_init = (X_init+X_init')/2;

X_init = X_init(:);

X_next = IBGP(X_init, affRnd, zeros(10,1), 0, 0.002, 400, 200, 2);
X_mat2 = reshape(X_next,[10 10]);
figure(2);
subplot(1,2,1), imshow(X_mat2,[]), title('result from seed 2');
subplot(1,2,2), imshow(reshape(X_init,[10 10]),[]), title('seed 2');
display(['affinity score IBGP seed 2: ',num2str(X_next'*affRnd*X_next)]);

%% RRWM method

% algpar.iterMax1 = 300;
% X_rrwm = RRWM(affRnd,10,10,algpar);
% X_cho = reshape(X_rrwm,[10 10]);
% % X_cho = (X_cho+X_cho')/2;
% 
% figure(3);
% subplot(1,3,1), imshow(X_mat1,[]), title('IBGP result');
% subplot(1,3,2), imshow(ones(10,10)*max(X_cho(:))-X_cho,[]), title('RRWM inverse result');
% subplot(1,3,3), imshow(X_cho,[]), title('RRWM result');

%% IPFP
[X,Y] = meshgrid(1:nodeCnt,1:nodeCnt);
X = X(:); Y=Y(:);
[sol, x_opt, score, score_sol]  = IPFP(affRnd, zeros(100,1), X_init_1, Y', X', 200);
display(['affinity score IPFP seed 1: ',num2str(x_opt'*affRnd*x_opt)]);
[sol, x_opt, score, score_sol]  = IPFP(affRnd, zeros(100,1), X_init_2(:), Y', X', 200);
display(['affinity score IPFP seed 2: ',num2str(x_opt'*affRnd*x_opt)]);
x_opt = reshape(x_opt,[10 10]);
figure(5);
subplot(1,2,1), imshow(X_mat2,[]), title('IBGP');
subplot(1,2,2), imshow(x_opt,[]), title('IPFP');