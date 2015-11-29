function LGPR_PLOT(X,Y,Interpolant,X_train,Z_prediction,S_prediction,ilgpr,include_centers,figureHandle,sp1,sp2,sp3)

    % each COLUMN of X_train is a training point    
    % Inperpolant is an interpolation object for the ground truth
    % X and Y are the complete grid of locations (created from a meshgrid function)
    % 

    addpath(genpath('./freezeColors'));
    addpath(genpath('./subplot_tight'));
    
    figure(figureHandle);
    
%     subplot(sp1,sp2,sp3);
    subplot_tight(sp1,sp2,sp3);
    
    Z_train = Interpolant(X_train(1,:),X_train(2,:));
    
    plot3(X_train(1,:),X_train(2,:),Z_train,'k+','MarkerSize',5,'LineWidth',3); % training data
    hold on
    freezeColors
    colormap('winter');
    surf(X,Y,reshape(Z_prediction,size(X)),'EdgeColor',[0.7 0.7 0.7],'FaceAlpha',0.75,'EdgeAlpha',0.2) % the prediction surface
    freezeColors
    colormap([0.8 0.8 0.8]);
    surf(X,Y,reshape(Z_prediction + 2*S_prediction,size(X)),'EdgeColor', [0.7 0.7 0.7],'FaceAlpha',0.5,'EdgeAlpha',0.2)
    surf(X,Y,reshape(Z_prediction - 2*S_prediction,size(X)),'EdgeColor', [0.7 0.7 0.7],'FaceAlpha',0.5,'EdgeAlpha',0.2)
    freezeColors
    axis([0 40 0 40 min(Z_train-10) max(Z_train+10)])
    
    if include_centers
        plot_ILGPR_centers(ilgpr,Interpolant,figureHandle);
    end
    
    hold off

end