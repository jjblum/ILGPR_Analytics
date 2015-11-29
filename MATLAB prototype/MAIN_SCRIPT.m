clear all
close all
clc

addpath(genpath('./gpml-v3.5')); startup;
addpath(genpath('./freezeColors'));

%% 1D sine wave example
% ground_truth_mean = 15;
% ground_truth = @(x) sin(x) + ground_truth_mean;
% noisy_func = @(x) ground_truth(x) + 0.5*randn(size(x));
% 
% x_all = (0:0.01:40)';
% x_fraction = 10;
% sparse_x = x_all(sort(randperm(length(x_all),round(length(x_all)/x_fraction)),'ascend'));
% y = noisy_func(sparse_x);
% 
% predictionX = x_all'; % note that each X is a column vector
% predictionZ = zeros(size(predictionX,1),1); % assume value is 0
% predictionS = 10*ones(size(predictionX,1),1); % assume +/- 10
% ilgpr = ILGPR(predictionX,predictionZ,predictionS); % the ILGPR object
% 
% N = length(sparse_x);
% datum = cell(size(predictionX,1),1); % training data
% for j = 1:N
%     x = sparse_x(j);
%     z = y(j);
%     datum{j} = Datum(x,z,j);
%     ilgpr.newDatum(datum{j});
% end
% 
% predictionZ = zeros(size(x_all));
% predictionS = zeros(size(x_all));
% for j = 1:length(x_all)
%     x = predictionX(j);
%     [predictionZ(j),predictionS(j)] = ilgpr.predict(x);
% end
% 
% temp = [ilgpr.LGPs{:}];
% temp = temp(vertcat(temp.started)==1);
% temp2 = [temp.hyp];
% hyp_cov = exp(horzcat(temp2.cov))
% PLOT_COUNT = 1 + length(temp);
% PLOT_COLS = 3;
% PLOT_ROWS = ceil(PLOT_COUNT/PLOT_COLS);
% 
% figure;
% subplot(PLOT_ROWS,PLOT_COLS,1)
% hold on
% % fill([predictionX fliplr(predictionX)],[predictionZ'+2*predictionS' fliplr(predictionZ'-2*predictionS')],'c')
% plot(sparse_x,y,'k+');
% plot(x_all,ground_truth(x_all),'b--');
% plot(predictionX,predictionZ,'r-','LineWidth',2);
% plot(predictionX,predictionZ + 2*predictionS,'c-','LineWidth',2)
% plot(predictionX,predictionZ - 2*predictionS,'c-','LineWidth',2)
% for j = 1:ilgpr.nLGPs
%     if ilgpr.LGPs{j}.started == 1
%         plot(ilgpr.LGPs{j}.u,ilgpr.LGPs{j}.predict(ilgpr.LGPs{j}.u),'ms','MarkerSize',20,'LineWidth',4);
%     end
% end
% hold off
% 
% n = 0;
% for j = 1:ilgpr.nLGPs
%     if ilgpr.LGPs{j}.started == 1
%         n = n+1;
%         subplot(PLOT_ROWS,PLOT_COLS,1+n)
%         LGP = ilgpr.LGPs{j};
%         
% %         if exp(LGP.hyp.cov(2)) < 0.1
% %             keyboard
% %         end
%         
%         
% %         figure
%         hold on
%         plot(x_all,ground_truth(x_all),'b--');
%         plot(LGP.X,LGP.Z+LGP.Zmean,'k+')
%         plot(LGP.u,LGP.predict(LGP.u),'ms','MarkerSize',20,'LineWidth',4);
%         local_x_all = linspace(min(LGP.X)-5,max(LGP.X)+5,100);
%         local_predictionZ = zeros(length(local_x_all),1);
%         local_predictionS = zeros(length(local_x_all),1);
%         for jj = 1:length(local_x_all)
%             [local_predictionZ(jj),local_predictionS(jj)] = LGP.predict(local_x_all(jj));
%         end
%         plot(local_x_all,local_predictionZ,'r-','LineWidth',2)
%         plot(local_x_all,local_predictionZ + 2*local_predictionS,'c-','LineWidth',2)
%         plot(local_x_all,local_predictionZ - 2*local_predictionS,'c-','LineWidth',2)
%         axis([min(x_all), max(x_all) -2+ground_truth_mean 2+ground_truth_mean])
%         title_string = sprintf('LGP # %d',j);
%         title(title_string,'FontSize',14);
%     end
% end
% 
% return;

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

%% 2D example, premade map
load('./map_0001')
heatmap_grid = heatmap_grid + 200; % offset upward to check algorithm mean prediction


N = 400; % maximum number of sample locations
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
    fprintf('Adding training datum # %d of %d\n',j,N); 
    x = Xz(:,j);
    z = myInterpolant(x(1),x(2));
    datum{j} = Datum(x,z,j);
    ilgpr.newDatum(datum{j});
end

temp = [ilgpr.LGPs{:}];
temp = temp(vertcat(temp.started)==1);
temp2 = [temp.hyp];
hyp_cov = exp(horzcat(temp2.cov))

[sMSE,Z_test,S_test] = LGPR_PREDICT(ilgpr,predictionX,myInterpolant);

fprintf('sMSE = %.4f\n',sMSE);

% w1 = figure;
% LGPR_PLOT(X,Y,myInterpolant,Xz,Z_test,S_test,ilgpr,1,w1,4,4,[1,2,5,6]);
% title('BEFORE Hyperparameter Re-Optimization','FontSize',20)

% re-optimize hyperparameters from scratch to make sure you don't have crazy GPs -- note that this may not be necessary
% for j = 1:ilgpr.nLGPs
%     ilgpr.LGPs{j}.reOptimizeHyperparameters();
% end
% w2 = figure;
% LGPR_PLOT(X,Y,myInterpolant,Xz,Z_test,S_test,ilgpr,1,w2);
% title('AFTER Hyperparameter Re-Optimization','FontSize',20)

% load junk

w3 = figure;
title(sprintf('sMSE = %.4f',sMSE));
LGPR_PLOT(X,Y,myInterpolant,Xz,Z_test,S_test,ilgpr,1,w3,4,4,[1,2,5,6]);


subplots = [3, 4, 7, 8, 9:16];

for j = 1:min(ilgpr.nLGPs,length(subplots))
    
    fprintf('Generating plot for LGP # %d of %d\n',j,min(ilgpr.nLGPs,length(subplots)));
    
    LGPs_temp = cell(1,1);
    LGPs_temp{1} = ilgpr.LGPs{j};
    Xz = ilgpr.LGPs{j}.getX();
    ilgpr_temp = ILGPR(predictionX,predictionZ,predictionS); % the ILGPR object
    ilgpr_temp.setLGPs(1,LGPs_temp);
    
    x_temp = LGPs_temp{1}.getX();
    bounding_box_edge = max( max(x_temp(1,:))-min(x_temp(1,:)) , max(x_temp(2,:))-min(x_temp(2,:)) );
    
    distances = pdist2(predictionX',LGPs_temp{1}.u');
    local_predictionX = predictionX(:,distances < bounding_box_edge/2);
    [localX,localY] = meshgrid(unique(local_predictionX(1,:)'),unique(local_predictionX(2,:)'));
    local_predictionX = horzcat(localX(:),localY(:))';
    
    [sMSE_temp,Z_test_temp,S_test_temp] = LGPR_PREDICT(ilgpr_temp,local_predictionX,myInterpolant);
    LGPR_PLOT(localX,localY,myInterpolant,Xz,Z_test_temp,S_test_temp,ilgpr_temp,1,w3,4,4,subplots(j));
end


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