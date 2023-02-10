function [checker]=sanity_check(L_cell,L_RVE,t_cell,t_matrix,t_matrix_layer,var_rnd_points,graph,n_points,n_points_per_line);
    checker=0;
if n_points <= 1 || n_points_per_line <= 1
    disp('Error: Number of points for voronoi must be greater than 1.')
    disp('Please change "L_cell" or "L_RVE".')
    disp('------------End of program------------')
    checker=1;
elseif t_matrix>L_RVE || t_matrix_layer>L_RVE
    disp('Error: Thickness of the matrix material can not be greater than the size of the RVE area.')
    disp('Please change "t_matrix", "t_matrix_layer" or "L_RVE".')
    disp('------------End of program------------')
    checker=1;
elseif t_matrix > L_cell*1.0
    disp('Error: Thickness of the matrix material should not be greater than 100% of the "L_cell".')
    disp('Please change at least one value.')
    disp('------------End of program------------')
    checker=1;
elseif var_rnd_points > 1 || var_rnd_points < 0
    disp('Error: The "var_rnd_points" has to be inside the interval [0 - 1].')
    disp('Please change the value.')
    disp('------------End of program------------')
    checker=1;
elseif n_points_per_line>10 || n_points > 100
    disp('!!!')
    disp('Hey, I guess, that this is totally not what you want! Stop this calculation now and try to reduce the points per line to maximum 20! Otherwise you will wait a lot of time.')
    disp('!!!')
end
if graph~=0 && graph~=1 && graph~=2
    disp('!!!')
    disp('Your input for the "graphic output" should be "0", "1" or "2". Now you will get no graphic output.')
    disp('!!!')
    graph=0;
end