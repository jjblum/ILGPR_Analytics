function [out, i_keep] = filter_area_pts(in_pts, ref_pts)
%% FILTER AREA POINTS
% Given a set of input points, filter out points that do not fall within
% the closed polygon described by ref_pts.
%
% This function is still developmental and has not been fully tested.

%%
if ref_pts(1,:) ~= ref_pts(end,:)
    ref_pts(end+1,:) = ref_pts(1,:);
end

%%
x_keep = zeros(1, length(in_pts));

min_x = min(ref_pts(:,1));
max_x = max(ref_pts(:,1));
min_y = min(ref_pts(:,2));
max_y = max(ref_pts(:,2));

%%
for i = 1:size(in_pts,1)
    xy = in_pts(i,:);
    
    if (xy(1) < min_x || xy(1) > max_x) && (xy(2) < min_y || xy(2) > max_y)
        continue
    end
    
    if any(xy(1) == ref_pts(:,1)) || any(xy(2) == ref_pts(:,2))
        continue
    end
    
    x_xchange = find(diff(xy(1) < ref_pts(:,1)));
    x_ychange = find(diff(xy(2) < ref_pts(:,2)));
    
    assert(mod(length(x_xchange),2) == 0);
    assert(mod(length(x_ychange),2) == 0);
    
    %%
    y_int_bounds = mean([ref_pts(x_xchange,2), ref_pts(x_xchange+1,2)],2);
    
    if (mod(sum(xy(2) > y_int_bounds),2) == 0) || ...
            (mod(sum(xy(2) < y_int_bounds),2) == 0)
        continue
    end
    
    x_int_bounds = mean([ref_pts(x_ychange,1), ref_pts(x_ychange+1,1)],2);
    
    if (mod(sum(xy(1) > x_int_bounds),2) == 0) || ...
            (mod(sum(xy(1) < x_int_bounds),2) == 0)
        continue
    end
    
    %%
    x_keep(i) = 1;
end

%%
i_keep = find(x_keep);
out = in_pts(i_keep,:);
