function LGPR_PLOT(X,Y,Interpolant,X_train,Z_prediction,S_prediction,ilgpr,include_centers,figureHandle)

    % each COLUMN of X_train is a training point    
    % Inperpolant is an interpolation object for the ground truth
    % X and Y are the complete grid of locations (created from a meshgrid function)
    % 

    addpath(genpath('./freezeColors'));
    
    figure(figureHandle);
    
    Z_train = Interpolant(X_train(1,:),X_train(2,:));
    
    plot3(X_train(1,:),X_train(2,:),Z_train,'k+','MarkerSize',5,'LineWidth',3); % training data
    hold on
    freezeColors
    colormap('copper');
    surf(X,Y,reshape(Z_prediction,size(X)),'EdgeColor',[0.7 0.7 0.7],'FaceAlpha',0.5) % the prediction surface
    freezeColors
    colormap([0.8 0.8 0.8]);
    surf(X,Y,reshape(Z_prediction + 2*S_prediction,size(X)),'EdgeColor', [0.7 0.7 0.7],'FaceAlpha',0.5)
    surf(X,Y,reshape(Z_prediction - 2*S_prediction,size(X)),'EdgeColor', [0.7 0.7 0.7],'FaceAlpha',0.5)
    freezeColors
    axis([-inf inf -inf inf min(Z_train-10) max(Z_train+10)])
    
    if include_centers
        plot_ILGPR_centers(ilgpr,Interpolant,figureHandle);
    end
    
    hold off

end