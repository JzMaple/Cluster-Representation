% config
n = 60;
sigma = 0.1;

% choose a grid of points
[x y] = meshgrid((1:n)/n);

% gaussian kernel
L = exp(- (bsxfun(@minus,x(:),x(:)').^2 + ...
           bsxfun(@minus,y(:),y(:)').^2) / sigma^2);

% sample
sample_dpp = sample_dpp(decompose_kernel(L));
sample_indep = randsample(n*n,length(sample_dpp));
  
% plot
subplot(1,2,1);
plot(x(sample_dpp),y(sample_dpp),'b.');
axis([0 1.02 0 1.02]);
axis square;
set(gca,'YTick',[]);
set(gca,'XTick',[]);
xlabel('DPP');

subplot(1,2,2);
plot(x(sample_indep),y(sample_indep),'r.');
axis([0 1.02 0 1.02]);
axis square;
set(gca,'YTick',[]);
set(gca,'XTick',[]);
xlabel('Independent');
