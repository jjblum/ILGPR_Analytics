clear all
close all
clc

addpath(genpath('./gpml-v3.5')); startup;

load('./map_0001')

% %% use GPML on a 2-D example, generated using GPML sampling
% ell = 10; sf = 15; sn = 0.1; % sn must be nonzero to avoid not-positive-definite errors in GPML code
% meanfunc = @meanConst; hyp.mean = 50;
% covfunc = @covSEiso; hyp.cov = log([ell; sf]);
% likfunc = @likGauss; hyp.lik = log(sn);
% 
% %% SAMPLE USING A PRIORI HYPERPARAMETERS
% n = 200; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NUMBER OF DATA POINTS
% x = 10*gpml_randn(rand(1), n, 2);
% K = feval(covfunc, hyp.cov, x);
% K = K + 0.00001*eye(size(K));
% [V,D] = eig(K);
% mu = feval(meanfunc, hyp.mean, x);
% y = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(rand(1), n, 1);
% 
% %% EXAMINE FUNCTION WITH A PRIORI HYPERPARAMETERS
% test_x1 = linspace(min(x(:,1))-1, max(x(:,1))+1, 101)';
% test_x2 = linspace(min(x(:,2))-1, max(x(:,2))+1, 101)';
% [X1,X2] = meshgrid(test_x1,test_x2);
% test_x = horzcat(X1(:),X2(:));
% [f_at_test_x, s2_at_test_x] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, test_x);
% F_at_test_x = reshape(f_at_test_x,size(X1));
% S2_at_test_x = reshape(s2_at_test_x,size(X1));
% F_plus_S = F_at_test_x + 2*sqrt(S2_at_test_x);
% F_minus_S = F_at_test_x - 2*sqrt(S2_at_test_x);
% 
% figure
% colormap copper
% shading interp
% hold on
% surf(X1,X2,F_at_test_x,'EdgeColor','none','FaceAlpha',1)
% surf(X1,X2,F_plus_S,'EdgeColor','none','FaceAlpha',0.25,'EdgeAlpha',0.25)
% surf(X1,X2,F_minus_S,'EdgeColor','none','FaceAlpha',0.25,'EdgeAlpha',0.25)
% plot3(x(:,1),x(:,2),y,'+','LineWidth',2)
% xlabel('x'); ylabel('y'); zlabel('z');
% title('A priori hyperparameters');
% set(gca, 'FontSize', 14)

%% APPLY ILGPR TO APPROXIMATE THE A PRIORI FULL GP
N = 100;
Xz = 10*gpml_randn(rand(1), N, 2);
ilgpr = ILGPR();
for j = 1:N
%     z = interp2(X,Y,heatmap_grid,min(X(:)) + rand(1)*(max(X(:))-min(X(:))),min(Y(:)) + rand(1)*(max(Y(:))-min(Y(:))));
    x = Xz(j,:)';
    a = LGP();

end