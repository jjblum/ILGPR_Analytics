function [sMSE,Z_test,S_test] = LGPR_PREDICT(ilgpr,X_test,Interpolant)

    % each COLUMN of X_test is a training point  

    N_test = size(X_test,2);
    test_error = zeros(size(X_test,2),1); % training error
    Z_test = zeros(size(X_test,2),1);
    S_test = zeros(size(X_test,2),1);
    for j = 1:N_test
%         fprintf('Evaluating test point # %d of %d\n',j,N_test);
        x = X_test(:,j);
        z = Interpolant(x(1),x(2));
        [Z_test(j),S_test(j)] = ilgpr.predict(x);
        test_error(j) = Z_test(j) - z;
    end
    true_z = Interpolant(X_test(1,:)',X_test(2,:)');
    sMSE = mean(test_error.^2)/std(true_z)^2;
    sMSE_baseline = mean((mean(true_z)-true_z).^2)/std(true_z)^2;

end