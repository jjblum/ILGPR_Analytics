function [varargout] = visualize_ilgpr(ilgpr, varargin)
%% VISUALIZE ILGPR
% visualize_ilgpr(ilgpr, x1, x2)
% predictions = visualize_ilgpr(ilgpr, x1, x2)
%
% Visualizes the output of Improved Local Gaussian Process Regression.
% 
% EXAMPLE USAGE:
% visualize_ilgpr(ilgpr, 1:10, 1:10);
%   - Plots the visualizations for the range of 1:10 in x and 1:10 in y.
% predictions = visualize_ilgpr(ilgpr, 1:10, 1:10);
%   - Plots the visualizations for the range of 1:10 in x and 1:10 in y,
%     and provides the prediction matrix as output.
%
% INPUTS:
%   - ilgpr: a trained ILGPR class
%   - x1: the data in the first dimension of x to visualize
%   - x2: the data in the second dimension of x to visualize
% OUTPUTS:
%   - None required
%
% OTHER NOTES:
%   - This function only visualizes 2-D data at this time.

%% Input argument handling
% Only allow 2 or 3 input arguments
narginchk(2,3);

% Check dimensionality
if length(varargin) == 1
    % 1-D Data
    D = 1;
    assert(0);  % not implemented yet
else
    % 2-D Data
    D = 2;
    
    % Unpack varargin
    x1 = varargin{1};
    x2 = varargin{2};
end

assert(isequal(size(x1, 1), size(x2, 1)), ...
    'x1 and x2 must be the same dimension size.');

% Constrain the visualization to work for 1-D or 2-D cases
assert(ismember(D, [1 2]), ...
    'This function can only visualize 1-D or 2-D cases at this time.');

%% Assemble grid of predictions
[X1,X2] = meshgrid(x1,x2);
prediction_grid = zeros(length(x1), length(x2));
for r = 1:size(prediction_grid,1)
    for c = 1:size(prediction_grid,2)
        x_sample = [x1(r); x2(c)];
        sample_datum = Datum(x_sample, [], []);
        prediction_grid(r,c) = ilgpr.predict(sample_datum);
    end
end

%% Get the location of the LGPs
loc_lgps = zeros(D, length(ilgpr.LGPs));
pred_lpgs = zeros(1, size(loc_lgps,2));
for i = 1:size(loc_lgps,2)
    loc_lgps(:,i) = ilgpr.LGPs{i}.u;
    sample_datum = Datum(loc_lgps(:,i), [], []);
    pred_lpgs(i) = ilgpr.predict(sample_datum);
end

%% Plot 3D prediction surface
f = figure; set(f, 'Position', [150 150 1280 800]);
h1 = surf(X1, X2, prediction_grid', 'EdgeColor', 'none');

axis image;

% Plot LGPs
hold on;
h2 = plot3(loc_lgps(1,:), loc_lgps(2,:), pred_lpgs, 'rx', 'MarkerSize', 8,...
    'LineWidth',2);

title('Prediction surface');
legend(h2,{'Local Gaussian Process'},'Location','Best');
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'FontSize',14);

%% Plot heat map contour
f = figure; set(f, 'Position', [150 150 1280 800]);
[C,h1] = contourf(X1, X2, prediction_grid');
clabel(C, h1, 'LineStyle', '--');

axis image;

% Plot LGPs
hold on;
h2 = plot(loc_lgps(1,:), loc_lgps(2,:), 'rx', 'MarkerSize', 8,...
    'LineWidth',2);

title('Prediction contour');
legend(h2,{'Local Gaussian Process'},'Location','Best');
xlabel('x');
ylabel('y');
set(gca,'FontSize',14);

%% Output argument handling
if nargout > 0
    varargout{1} = prediction_grid;
end
