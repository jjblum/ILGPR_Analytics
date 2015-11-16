clear all
close all
clc

addpath(genpath('./gpml-v3.5')); startup;

load('./map_0001')

%% use GPML on a 2-D example, generated using GPML sampling
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
N = 400; % maximum number of sample locations
Xz = 20+5*gpml_randn(rand(1), N, 2)'; % predetermined random set of sample locations, note that each X is a column vector
Xz(Xz > 40) = 40;
Xz(Xz < 1) = 1;
predictionX = heatmap(:,1:2)'; % note that each X is a column vector
predictionZ = zeros(size(predictionX,1),1); % assume value is 0
predictionS = 10*ones(size(predictionX,1),1); % assume +/- 10
ilgpr = ILGPR(predictionX,predictionZ,predictionS); % the ILGPR object
myInterpolant = griddedInterpolant(X,Y,heatmap_grid,'cubic');

datum = cell(size(predictionX,1),1); % training data
for j = 1:N
    x = Xz(:,j);
    z = myInterpolant(x(1),x(2));
    datum{j} = Datum(x,z,j);
    ilgpr.newDatum(datum{j});
end

%% Test prediction on training set
train_label = zeros(N,1); % trainting label
train_pred = zeros(N,1); % trainting prediction
for j = 1:N
    train_label(j) = datum{j}.getZ();
    train_pred(j) = ilgpr.predict(datum{j});
end
train_error = train_label-train_pred; % training error
train_mse = mean((train_label-train_pred).^2); % MSE of training data
train_nmse = train_mse/mean(train_label)/mean(train_pred); % nMSE of trainting data
display(train_nmse);

%% Test prediction on testing set
N_test = 200; % Number of testing data
Xz_test = 20+5*gpml_randn(rand(1), N_test, 2)';
Xz_test(Xz_test > 40) = 40;
Xz_test(Xz_test < 1) = 1;

test_label = zeros(N_test,1); % testing label
test_pred = zeros(N_test,1); %testing prediction
for j = 1:N_test
    x = Xz_test(:,j);
    test_label(j) = myInterpolant(x(1),x(2));
    datum = Datum(x,z,j);
    test_pred(j) = ilgpr.predict(datum);
end
test_error = test_label-test_pred; % testing error
test_mse = mean((test_label-test_pred).^2); % MSE of testing data
test_nmse = test_mse/mean(test_label)/mean(test_pred); % nMSE of testing data
display(test_nmse);

%% visualization
% temp = [ilgpr.LGPs{:}];
% loc_lgps = [temp.u]; clear temp;
% 
% res = 100;
% x_1 = linspace(min([Xz_test(1,:),loc_lgps(1,:)]), ...
%                max([Xz_test(1,:),loc_lgps(1,:)]), res+1);
% x_2 = linspace(min([Xz_test(2,:),loc_lgps(2,:)]), ...
%                max([Xz_test(2,:),loc_lgps(2,:)]), res+1);
% predictions = visualize_ilgpr(ilgpr, x_1, x_2);