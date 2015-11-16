clear all
close all
clc

addpath(genpath('./gpml-v3.5')); startup;
addpath(genpath('./freezeColors'));

load('./map_0002')

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
N = 500; % maximum number of sample locations
Xz = 20+10*gpml_randn(rand(1), N, 2)'; % predetermined random set of sample locations, note that each X is a column vector
Xz(Xz > 40) = 40;
Xz(Xz < 1) = 1;
predictionX = heatmap(:,1:2)'; % note that each X is a column vector
predictionZ = zeros(size(predictionX,1),1); % assume value is 0
predictionS = 10*ones(size(predictionX,1),1); % assume +/- 10
ilgpr = ILGPR(predictionX,predictionZ,predictionS); % the ILGPR object
myInterpolant = griddedInterpolant(X,Y,heatmap_grid,'cubic');

datum = cell(size(predictionX,1),1); % training data
train_error = zeros(size(predictionX,1),1); % training error
for j = 1:N
    x = Xz(:,j);
    z = myInterpolant(x(1),x(2));
    datum{j} = Datum(x,z,j);
    ilgpr.newDatum(datum{j});
end

% training set prediction
predictionZ = zeros(N,1);
predictionS = zeros(N,1);
for j = 1:N    
    x = datum{j}.getX();
    [predictionZ(j),predictionS(j)] = ilgpr.predict(x);
    train_error(j) = datum{j}.getZ() - predictionZ(j);
%     disp(sprintf('Training Evaluation, Data at [%f,%f] = %f, Prediction = %f, Error = %f\n',x(1),x(2),datum{j}.getZ(),predictionZ(j),train_error(j)));
end
surf(X,Y,heatmap_grid); hold on;
plot3(Xz(1,:),Xz(2,:),predictionZ,'k+','MarkerSize',5,'LineWidth',3);
% plot3(Xz(1,:),Xz(2,:),predictionS,'+','Color',[0.7 0.7 0.7],'MarkerSize',5,'LineWidth',3);
for j = 1:ilgpr.nLGPs
    plot3(ilgpr.LGPs{j}.u(1),ilgpr.LGPs{j}.u(2),myInterpolant(ilgpr.LGPs{j}.u(1),ilgpr.LGPs{j}.u(2)),'ms','MarkerSize',10,'LineWidth',5);
end

% test prediction on entire grid
Xz_test = predictionX;
N_test = size(Xz_test,2);
test_error = zeros(size(Xz_test,2),1); % training error
test_Z = zeros(size(Xz_test,2),1);
test_S = zeros(size(Xz_test,2),1);
for j = 1:N_test
    x = Xz_test(:,j);
    z = myInterpolant(x(1),x(2));
    datum_test = Datum(x,z,j);
    [test_Z(j),test_S(j)] = ilgpr.predict(datum_test.getX());
    test_error(j) = test_Z(j) - z;
end

figure
plot3(Xz(1,:),Xz(2,:),myInterpolant(Xz(1,:),Xz(2,:)),'k+','MarkerSize',5,'LineWidth',3); % training data
hold on
freezeColors
colormap('winter');
surf(X,Y,reshape(test_Z,size(X)),'EdgeColor',[0.7 0.7 0.7]) % the prediction surface
freezeColors
colormap([0.8 0.8 0.8]);
surf(X,Y,reshape(test_Z + 2*test_S,size(X)),'EdgeColor', [0.7 0.7 0.7],'FaceAlpha',0.5)
surf(X,Y,reshape(test_Z - 2*test_S,size(X)),'EdgeColor', [0.7 0.7 0.7],'FaceAlpha',0.5)
freezeColors
hold off


% figure
% hold on
% for j = 1:N_test
%     error_ratio = abs(test_error(j))/max(abs(test_error(:)));
%     if sign(test_error(j)) < 0
%         color = [1-error_ratio 1-error_ratio 1];
%     else
%         color = [1 1-error_ratio 1-error_ratio];
%     end
%     plot3(Xz_test(1,j),Xz_test(2,j),test_error(j),'.','MarkerSize',20,'Color',color);
% end
% hold off
% 
% % visualization
% temp = [ilgpr.LGPs{:}];
% loc_lgps = [temp.u]; clear temp;
% 
% res = 20;
% x_1 = linspace(min([Xz_test(1,:),loc_lgps(1,:)]), ...
%                max([Xz_test(1,:),loc_lgps(1,:)]), res+1);
% x_2 = linspace(min([Xz_test(2,:),loc_lgps(2,:)]), ...
%                max([Xz_test(2,:),loc_lgps(2,:)]), res+1);
% predictionZ = visualize_ilgpr(ilgpr, x_1, x_2);