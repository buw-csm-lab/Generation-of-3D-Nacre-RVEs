function [list_of_rnd_points, periodic_point_list]=RasterPointGeneration(L_RVE,n_points,L_cell,var_rnd_points);

disp('Raster point generation...');


list_of_rnd_points=[];

% Create an even raster of points
for m=1:1:n_points
    x=(m/n_points)*L_RVE;
    for n=1:1:n_points
        y=(n/n_points)*L_RVE;
         if rem(m,2)    %give every secound column a
             y=y-(1/(2*n_points))*L_RVE;    %offset of -1/2*L_RVE
         end
        point=[x,y];
        %append the currently generated point to list_of_rnd_points
        list_of_rnd_points=[list_of_rnd_points;point];
    end
end

% Scale the variation of position of random points with L_cell
var_rnd_points=var_rnd_points*L_cell;

% Each coordinate of the list_of_rnd_points gets decreased by a random
% value multiplied with the parameter 'var_rnd_points'
for i=1:1:length(list_of_rnd_points)
    %multiply the maximal displacement of voronoi-seeds with a random number between 0 and 1
    x_var_rnd_points=(var_rnd_points).*rand;
    y_var_rnd_points=(var_rnd_points).*rand;
    %subtract the displacement from the current position of points 
    list_of_rnd_points(i,1)=list_of_rnd_points(i,1)-x_var_rnd_points;
    list_of_rnd_points(i,2)=list_of_rnd_points(i,2)-y_var_rnd_points;
end

%create 'periodic points' from 'list_of_rnd_points' in the 8 surrounding areas
[periodic_point_list,number_of_points]=PeriodicPoints(L_RVE,list_of_rnd_points);
fprintf('%d points were placed in the RVE area.\n',number_of_points)