function [list_of_rnd_points, periodic_point_list]=RSAPointGeneration(L_RVE,min_d,abort_criterion_RSA);

disp('RSA Point generation...');

% Parameter for the stand-alone version
% L_RVE=13.;
% min_d=2;
% critical_value=100000;
% A_cell=sqrt(3)*3/2.*(min_d/2.)^2;
% A_RVE=L_RVE^2;
% n_points=round(A_RVE/A_cell)

%create a random start point within the RVE and add it to the
%'list_of_rnd_points'
start_rnd_point=rand(1,2)*L_RVE;
list_of_rnd_points=[start_rnd_point];

%create 'periodic points' from 'list_of_rnd_points' in the 8 surrounding areas
[periodic_point_list,number_of_points]=PeriodicPoints(L_RVE,list_of_rnd_points);

securety_counter=0;
ii=1;

% while number of tries to place a new point is smaller than the abort_criterion
while securety_counter<abort_criterion_RSA

    %create a random point within the RVE'
    rnd_point=rand(1,2)*L_RVE;
    list_of_rnd_points_periodic_inc=[rnd_point;periodic_point_list];
    
    %calculate the distances between all fixed points and the new point
    list_of_distances = pdist(list_of_rnd_points_periodic_inc);
    size_of_d=size(list_of_distances,2);

    i=1;
    % for every distance between two points...
    while i<=size_of_d
        % if any distance is smaller than the minimum-distance-criterium, while loop
        % ends and a new point will be generated
        if list_of_distances(1,i) < min_d
           i=size_of_d+1;
           securety_counter=securety_counter+1;
        else
            % after comparing the minimum-distance-criterium with every distance...
            if i==size_of_d
                % ...the point generation was successfull
                % add the new point to the list of already fixed points('list_of_rnd_points')
                list_of_rnd_points= [rnd_point;list_of_rnd_points];
                
                %create 'periodic points' from 'list_of_rnd_points' in the 8 surrounding areas
                [periodic_point_list,number_of_points]=PeriodicPoints(L_RVE,list_of_rnd_points);
                
                securety_counter=0;
            end
            i=i+1;
        end
        % if the value of tries to place a new random point is higher than
        % the abort-criterion, then stop the function
        if securety_counter>=abort_criterion_RSA
            break
        end
    end
end
fprintf('%d points were placed in the RVE area.\n',number_of_points)
