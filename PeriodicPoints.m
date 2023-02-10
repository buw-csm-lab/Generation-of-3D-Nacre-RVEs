function [periodic_point_list,number_of_points]=PeriodicPoints(length_original_area,point_list);

% Copy every point within 'point_list' into the 8 areas, which are
% surrounding the quadratic area with the edge length
% 'length_original_area'
periodic_point_list=[];
number_of_points=length(point_list(:,1));
for i = -1:1:1
    for j = -1:1:1
        for k = 1:1:number_of_points
                a=point_list(k,1)+i*length_original_area;
                b=point_list(k,2)+j*length_original_area;
                periodic_point_list=[periodic_point_list;[a,b]];
        end
    end
end