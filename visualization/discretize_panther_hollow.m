%% DISCRETIZE PANTHER HOLLOW LAKE
% This script discretizes Panther Hollow Lake into a grid of GPS
% coordinate points. It saves the grid points in the file 
% "panther_hollow_gps_points.txt".

%%
clear all; close all; clc;

%% 
fid = fopen('panther_hollow_gps_boundary.txt');
assert(fid >= 3);
c_data = textscan(fid,'%f,%f');
fclose(fid);

gps = cell2mat(c_data);
longitude = gps(:,2);
latitude = gps(:,1);

%%
min_lat = min(latitude);
max_lat = max(latitude);
min_lon = min(longitude);
max_lon = max(longitude);

del = 5.04e-5;

lat_grid = min_lat:del:max_lat;
lon_grid = min_lon:del:max_lon;

grid_pts_all = [];
%
for x = 1:length(lon_grid)
    for y = 1:length(lat_grid)
        grid_pts_all(end+1,:) = [lon_grid(x), lat_grid(y)];
    end
end

boundary = [longitude, latitude];
grid_pts = filter_area_pts(grid_pts_all, boundary);

f = figure; set(f,'Position',[300 300 1280 800]);
h1 = plot(longitude, latitude, '.-');
axis equal; hold on; grid on;
h2 = plot(grid_pts(:,1), grid_pts(:,2), 'x');
xlabel('Longitude');
ylabel('Latitude');

title('Panther Hollow Lake - GPS Specifications');
set(gca,'FontSize',14);

legend([h1 h2],{'Boundary','Discretized interior'},'Location','SouthEast');

if 0
    saveas(f,'panther_hollow_gps','fig');
    saveas(f,'panther_hollow_gps','png');
end

%%
fid = fopen('panther_hollow_gps_points.txt','w');
A = [grid_pts(:,2), grid_pts(:,1)]';
fprintf(fid,'%6.8f,%6.8f\r\n',A);
fclose(fid);
