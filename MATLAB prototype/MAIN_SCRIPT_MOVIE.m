clear all
close all
clc

addpath(genpath('./gpml-v3.5')); startup;
addpath(genpath('./freezeColors'));

%% Load boat query points
fid = fopen('boat_query_points.txt');
assert(fid >= 3);
c_data = textscan(fid,'%f,%f');
fclose(fid);

boat_gps = cell2mat(c_data);
boat_gps = [boat_gps(:,2), boat_gps(:,1)];

% New code to increase the number of boat points
temp = zeros(size(boat_gps,1) + size(boat_gps,1) -1, 2);
temp(1:2:end,:) = boat_gps;
temp(2:2:end,:) = boat_gps(1:(end-1),:) + diff(boat_gps)/2;
boat_gps = temp;

% Determine max/min lat/long
del = 0.0025;
min_lon = -79.949355;
max_lon = min_lon + del;
max_lat = 40.437488;
min_lat = max_lat - del/2;

m_lat = (max_lat - min_lat)/(40 - 1);
b_lat = max_lat - m_lat*40;

m_lon = (max_lon - min_lon)/(40 - 1);
b_lon = max_lon - m_lon*40;

boat_points = [(boat_gps(:,1) - b_lon)./m_lon, ...
               (boat_gps(:,2) - b_lat)./m_lat];

%% Load all grid query points
fid = fopen('panther_hollow_gps_points.txt');
assert(fid >= 3);
c_data = textscan(fid,'%f,%f');
fclose(fid);

grid_pts = cell2mat(c_data);
grid_pts = [grid_pts(:,2), grid_pts(:,1)];

query_points = [(grid_pts(:,1) - b_lon)./m_lon, ...
                (grid_pts(:,2) - b_lat)./m_lat];

%% 2D example, premade map
load('./map_0001')
heatmap_grid = heatmap_grid + 200; % offset upward to check algorithm mean prediction


% N = 60; % maximum number of sample locations
% Xz = 20+10*gpml_randn(rand(1), N, 2)'; % predetermined random set of sample locations, note that each X is a column vector
% Xz(Xz > 40) = 40;
% Xz(Xz < 1) = 1;

Xz = boat_points';
N = size(Xz,2);

predictionX = heatmap(:,1:2)'; % note that each X is a column vector
predictionZ = zeros(size(predictionX,1),1); % assume value is 0
predictionS = 10*ones(size(predictionX,1),1); % assume +/- 10
ilgpr = ILGPR(predictionX,predictionZ,predictionS); % the ILGPR object
myInterpolant = griddedInterpolant(X,Y,heatmap_grid,'cubic');

datum = cell(size(predictionX,1),1); % training data
train_error = zeros(size(predictionX,1),1); % training error

sMSE = zeros(N,1);
Z_test = cell(N,1);
S_test = cell(N,1);

ALL_WEIGHTS = zeros(size(query_points,1),N);
ALL_S_VARS = zeros(size(ALL_WEIGHTS));

F(N) = struct('cdata',[],'colormap',[]); % make a movie
% video_writer = VideoWriter('./2D_example_movie.avi','Uncompressed AVI');
% open(video_writer);

nLGPs_all = zeros(1,N);

for j = 1:N
    fprintf('Adding training datum # %d of %d\n',j,N); 
    x = Xz(:,j);
    z = myInterpolant(x(1),x(2));
    datum{j} = Datum(x,z,j);
    ilgpr.newDatum(datum{j});
    
    [sMSE(j),Z_test{j},S_test{j}] = LGPR_PREDICT(ilgpr,predictionX,myInterpolant);
    
    w3 = figure;
%     set(w3,'Visible','off');
    set(w3,'Position',[70 1 1600 1200]);
    set(w3,'color',[1 1 1]);
%     mTextBox = uicontrol('style','text');
%     set(mTextBox,'String',sprintf('sMSE = %.4f',sMSE(j)));
%     set(mTextBox,'FontSize',30);
%     set(mTextBox,'Position',[1200 950 300 50]);
%     set(mTextBox,'BackgroundColor',[.9 .9 .9]);

    LGPR_PLOT(X,Y,myInterpolant,Xz(:,1:j),Z_test{j},S_test{j},ilgpr,1,w3,4,4,[1,2,5,6]);
    freezeColors;        

    subplots = 7:16;

    % sort the LGPs by the number of points assigned to them
    temp = [ilgpr.LGPs{:}];
    all_N = [temp.N];
    [sorted_N, sorted_N_index] = sort(all_N,'descend');

    for g = 1:min(ilgpr.nLGPs,length(subplots))

%         fprintf('Generating plot for LGP # %d of %d\n',j,min(ilgpr.nLGPs,length(subplots)));

        LGPs_temp = cell(1,1);
        LGPs_temp{1} = ilgpr.LGPs{sorted_N_index(g)};
                        
        X_train = ilgpr.LGPs{sorted_N_index(g)}.getX();
        ilgpr_temp = ILGPR(predictionX,predictionZ,predictionS); % the ILGPR object
        ilgpr_temp.setLGPs(1,LGPs_temp);

        x_temp = LGPs_temp{1}.getX();
        bounding_box_edge = max( max(x_temp(1,:))-min(x_temp(1,:)) , max(x_temp(2,:))-min(x_temp(2,:)) );

        distances = pdist2(predictionX',LGPs_temp{1}.u');
        local_predictionX = predictionX(:,distances < bounding_box_edge/2);
        [localX,localY] = meshgrid(unique(local_predictionX(1,:)'),unique(local_predictionX(2,:)'));
        local_predictionX = horzcat(localX(:),localY(:))';

        [~,Z_test_temp,S_test_temp] = LGPR_PREDICT(ilgpr_temp,local_predictionX,myInterpolant);
        LGPR_PLOT(localX,localY,myInterpolant,X_train,Z_test_temp,S_test_temp,ilgpr_temp,1,w3,4,4,subplots(g));
    end    
        
    subplot_tight(4,4,[3 4])
    plot(1:j,log(sMSE(1:j)),'rx-','LineWidth',3,'MarkerSize',5);
    hold on
    plot([1 j+1],[0 0],'b--','LineWidth',3);
    hold off
    set(gca,'FontSize',14)    
    axis([1 j+1 -Inf max(log(sMSE(1:j)))]);
    
    
    F(j) = getframe(w3); % capture a frame for the movie
%     writeVideo(video_writer,F(j));
    
    close all
    
    nLGPs_all(j) = ilgpr.nLGPs;
    
    % Query points
% %     weights = zeros(size(query_points,1),1);
% %     s_vars = zeros(size(weights));
% %     for q = 1:size(query_points,1)
% %         x = query_points(q,1);
% %         y = query_points(q,2);
% %         
% %         % find nearest 
% %         x_sweep = X(1:end,1);
% %         y_sweep = Y(1,1:end);
% %         
% %         [~,x_index] = min(abs(x_sweep - x));
% %         [~,y_index] = min(abs(y_sweep - y));
% %         
% %         ZMAT = reshape(Z_test{j},size(X));
% %         SMAT = reshape(S_test{j},size(X));
% %         
% %         % save weight
% %         weights(q) = ZMAT(x_index,y_index);
% %         s_vars(q) = SMAT(x_index, y_index);
% %     end
    [~,~,s_vars] = LGPR_PREDICT(ilgpr,query_points',myInterpolant);
    
    %ALL_WEIGHTS(:,j) = weights;
    ALL_S_VARS(:,j) = s_vars;
    
end

% close(video_writer);

% movie(F,1); % play the movie

fprintf('sMSE = %.4f\n',sMSE);

%% Dump to file
% ALL_S_VARS_reweighted = ALL_S_VARS./max(max(ALL_S_VARS))*3;
for j = 1:N
    fid = fopen(['uncertainty_',num2str(j),'.txt'],'w');
    weightdata = ALL_S_VARS(:,j)/max(S_test{j})*1;
    %A = [grid_pts(:,2), grid_pts(:,1), ALL_S_VARS(:,j)]';
    A = [grid_pts(:,2), grid_pts(:,1), weightdata]';
    fprintf(fid,'%6.8f,%6.8f,%1.3f,\r\n',A);
    fclose(fid);
    
    fid = fopen(['boat_',num2str(j),'.txt'],'w');
    boat_coords = boat_gps(j,:);
    fprintf(fid,'%6.8f,%6.8f',fliplr(boat_coords));
    fclose(fid);
end