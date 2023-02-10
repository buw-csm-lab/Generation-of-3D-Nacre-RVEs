%% GENERATION OF NACRE RVE

% Generates layers of nacre based on Voronoi tesselation with
% matrix of t_matrix thickness between the platelets.
% 
% number of generated layers: n_layer = (L_RVE/(t_cell+t_matrix_layer) 
%
% External requirements:
%   - Matlab AddOn "Mapping Toolbox"
%   - inpoly2.m | INPOLY: A fast points-in-polygon test by Darren Engwirda
%     tested last with Version 3.0.0.0 (https://github.com/dengwirda/inpoly)
%   - VoronoiLimit.m | VoronoiLimit(vararg​in) by Jakob Sievers 
%     tested last with Version 3.0.2.2 (https://de.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit-varargin)
%   - InterX.m | Curve intersections by NS
%     tested last with Version 1.5.0.0 (https://de.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections)
%
% Internal requirements / functions defined in files:
%       Cut_and_Close_a_Cell.m, PeriodicPoints.m, 
%		RasterPointGeneration.m, RSAPointGeneration.m, sanity_check.m
%
% output is saved in ./_data folder
% after generation run "a02_Nacre_FULL_3D_RVE.py" in ABAQUS in the same
% working directory.

disp('------------Start of program------------')
clear

%% Input section

% Define the user interface
prompt = {'L_cell:                  length of ceramic cells',...
    'L_RVE:                 length of RVE borders',...
    't_cell:                    thickness of ceramic cells',...
    't_matrix:                thickness of matrix between the cells in a layer',...
    't_matrix_layer:      thickness of matrix between the layers',...
    'method_of_pointgeneration: 0 = RSA, 1= Raster',...
    'var_rnd_points [0-1] (only necessary for Raster point generation)',...
    'graphic output:   no output = 0; Detail plot = 1;  2D RVE plot = 2'};
UI_title = 'User parameters';
dims = [1 70];
definput = {'2.5','10','0.9','0.2','0.2','0','0.45','2'};
answer = inputdlg(prompt,UI_title,dims,definput);

% Average cell size in [µm]
% reality: 5-8 µm
L_cell=2.3;
L_cell=str2num(strrep(answer{1,1},',','.'));

% Length of the RVE in [µm]
L_RVE=10;
L_RVE=str2num(strrep(answer{2,1},',','.'));

% Layer thickness of cells
% reality: 0.2-0.9 µm
t_cell=2;
t_cell=str2num(strrep(answer{3,1},',','.'));

% Thickness of matrix material in [µm]
% reality: ~0.02 µm
t_matrix=0.02;
t_matrix=str2num(strrep(answer{4,1},',','.'));

% Thickness of matrix material between each layer in [µm]
% reality: ~0.02 µm
t_matrix_layer=0.02;
t_matrix_layer=str2num(strrep(answer{5,1},',','.'));

% Decide the funktion of point generation
% 0 = RSA point generation
% 1 = Raster point generation
method_of_pointgeneration=0;
method_of_pointgeneration=str2num(strrep(answer{6,1},',','.'));

% Factor to displace the rnd-points before building cells with voronoi. (1=
% 100%=displacement of value from average cell size.) [standard=0.35]
var_rnd_points=0.35;
var_rnd_points=str2num(strrep(answer{7,1},',','.'));

% Different options for graphical output
% 0=no output 1=detail graphic; 2=2D-RVE;
graph=0;
graph=str2num(strrep(answer{8,1},',','.'));

if L_cell <= 0 || t_matrix <= 0 || L_RVE <= 0 || t_matrix_layer <= 0
    disp('Error: Input parameters "L_cell", "t_matrix", "n_points", "L_RVE", "t_matrix_layer" must be greater than zero.')
    disp('Please change those incorrect value(s).')
    return
end

%% Core data

% Calculate the number of layers
n_layer=(L_RVE/(t_cell+t_matrix_layer));
n_layer=round(n_layer);

% Set the borders of the 2D RVE area for plotting the RVE in MatLab
area=[0 0;0 L_RVE; L_RVE L_RVE;L_RVE 0;0 0];

% Scale axes for different output diagrams in MatLab
scale_low =      -0.1*L_RVE;
scale_low_RVE =  -0.1*L_RVE;
scale_low =      -L_RVE*1.1;
scale_high =      1.1*L_RVE;
scale_high_RVE =  1.1*L_RVE;
scale_high =      2.1*L_RVE;

% global parameter
number_of_element=0;
number_of_old_element=0;

% global counter
ccc=0;
cc=1;

% Some output values
fprintf('Number of layers - n_layer: %d \n\n',n_layer)
fprintf('Thickness of matrix material - t_matrix: %f µm\n\n',t_matrix)
fprintf('Thickness of matrix material between each layer - t_matrix_layer: %f µm\n\n',t_matrix_layer)
fprintf('Thickness of cell - t_cell: %f µm\n\n',t_cell)

% Calculate the true t_cell for building a cube
t_cell=(L_RVE/n_layer)-t_matrix_layer;
fprintf('Adapted thickness of cell - t_cell: %f µm\n\n',t_cell)
fprintf('Average size of cells: %f µm\n\n',L_cell)

% Calculate the number of points per line, which will be used to build
% polygons with voronoi
if method_of_pointgeneration==0
    % Calculate the value of of a ideal 
    A_cell=sqrt(3)*3/2.*(L_cell/2.)^2;
    A_RVE=L_RVE^2;
    n_points=round(A_RVE/A_cell);    
    
    % Assign arbitrary value for initialisation
    n_points_per_line=2;
    
elseif method_of_pointgeneration==1
    n_points_per_line=(L_RVE/L_cell);

    % Calculate the number of points per line
    % needs to be even integer
    n_points_r=floor(n_points_per_line);    %rounding off n_points
    if rem(n_points_r,2)    %if n_points_r is an uneven number... 
        n_points_per_line=ceil(n_points_per_line);  %...then round up n_points
    else    %if n_points_r is an even number...
        n_points_per_line=n_points_r;   %...use n_points_r
    end
    
    % if n_points_r is an uneven integer number (like 1.0, 3.0,..)
    if rem(n_points_per_line,2) %if n_points is still uneven...
        n_points_per_line=n_points_per_line+1;  %...add 1 to make it even
    end
    
    % Calculate the true average distance between two points in a line
    L_cell=L_RVE/n_points_per_line;
    fprintf('Number of points per line: %f \n\n',n_points_per_line)
    
    % Assign arbitrary value for initialisation
    n_points=16;
end

fprintf('Adapted average size of cells: %f µm\n\n',L_cell)


%% Sanity check [optional]

% Call the function "sanity_check" to check for correctness of the input
% parameters
[checker]=sanity_check(L_cell,L_RVE,t_cell,t_matrix,t_matrix_layer,var_rnd_points,graph,n_points,n_points_per_line);
if checker == 1
   return 
end

%% Main program
%Runs the code for each layer which will be built

for num_lay = 1:1:n_layer
    fprintf('\n*** Start of Layer %d / %d *** \n',num_lay,n_layer)
    
    % Set some output values for some plots
    if graph==1
        f1=figure('Name','Full finished layer');
        set(f1,'Visible','off')
    elseif graph==2
        figure('Name','2D RVE')
    elseif graph==0
        ff=figure('Name','graph==0','Visible','off');
    end

    %% RSA - Point generation
    if method_of_pointgeneration == 0
        
        %Set the abort criterion for tries to find a correct location for a
        %point
        abort_criterion_RSA=70000;
        
        % Call function for point generation method by RSA
        [list_of_rnd_points,periodic_point_list]=RSAPointGeneration(L_RVE,L_cell,abort_criterion_RSA);

        %% Raster - Point generation
    elseif method_of_pointgeneration==1
 
        % Call function for point generation method in a raster
        [list_of_rnd_points,periodic_point_list]=RasterPointGeneration(L_RVE,n_points_per_line,L_cell,var_rnd_points); 
    end
    
    
    %% Get geometrical information from voronoid-cmd
    
    % Define the area for Voronoi tesselation
    area_big=[-2*L_RVE,-2*L_RVE;-2*L_RVE,3*L_RVE;3*L_RVE,3*L_RVE;3*L_RVE,-2*L_RVE;-2*L_RVE,-2*L_RVE];
    
    % Provides the Voronoi tesselation of periodic_point_list 
    % V contains all vertices
    % C is a matrix of points used to build up each polygon
    % XY is not needed.
    [V,C,XY]=VoronoiLimit(periodic_point_list(:,1),periodic_point_list(:,2),'figure','off','bs_ext',area_big);
    for i=1:1:length(C)
        C{i}=(C{i})';
    end

    % Display C and V in a readable way 
%     for i = 1:length(C)
%         disp(C{i});
%     end
%     disp(V)

    %% Build up 2D-RVE (single wall cells)
    
    % This list contains in each row the position information of 
    % points used for one polygon
    list_of_polygons={};
    
     % A list for all vertices of all polygons
    list_of_elements=[];    
    
    list_of_old_polygons={};
    list_of_old_elements=[];
    counter_list_of_old_elements=1;
    list_old_to_new_vertex=[];
    
    
    disp('Build up 2D RVE with inner polygons...')
    % This loop runs for every generated polygon separated
    for i = 1:length(C)
        
        % A list for all vertices of the current polygon
        list_vertices_old_polygon=[];
        
        % A list for all elements=(edges of current polygon || lines between two
        % vertices of current polygon)
        list_of_current_elements=[];
        list_of_old_current_elements=[];

        
        % C{i} contains the position of the information about the
        % coordinates of each vertex of the current polygon
        for j=1:length(C{i})
            c=C{i};
            vx=V(c(j),1);
            vy=V(c(j),2);
            
            %added for output of one_walled cells
            if j~=length(C{i})
                vx_1=V(c(j+1),1);
                vy_1=V(c(j+1),2);
            else
                vx_1=V(c(1),1);
                vy_1=V(c(1),2);
            end
            % Write each x- and y-coordinate of every vertex of the
            % current polygon into the list_vertices_old_polygon
            list_vertices_old_polygon=[list_vertices_old_polygon;vx vy];
            
%             %added for output of one_walled cells
%             list_of_old_elements=[list_of_old_elements;counter_list_of_old_elements,vx,vy,vx_1,vy_1];
%             counter_list_of_old_elements=counter_list_of_old_elements+1;
            
            % Add the first vertex at the last position of the 
            % list_vertices_old_polygon. Necessary to line up the full
            % polygon later.
             if length(C{i})==j
                vx=V(c(1),1);
                vy=V(c(1),2);
                list_vertices_old_polygon=[list_vertices_old_polygon;vx vy];
            end
        end
        
%         for j=1:length(list_vertices_old_polygon)-1
%             list_of_old_elements=[list_of_old_elements;counter_list_of_old_elements,list_vertices_old_polygon(j,1),list_vertices_old_polygon(j,2),list_vertices_old_polygon(j+1,1),list_vertices_old_polygon(j+1,2)];
%             counter_list_of_old_elements=counter_list_of_old_elements+1;
%         end



        
        %disp('Create, cut and close single-wall polygons...')
        %% Call function for cutting and closing each inner cell at RVE border
        intersection_with_neighbour=1;
        [list_of_old_current_elements,number_of_old_element, list_of_old_current_polygons,list_of_old_elements]=...
            Cut_and_Close_a_Cell(list_vertices_old_polygon, list_of_old_current_elements,...
            number_of_old_element, area, L_RVE,list_of_old_elements,intersection_with_neighbour);
        intersection_with_neighbour=0;
        
        % Handle the list_of_current_polygons to the global list_of_polygons
        if list_of_old_current_polygons ~= 0
            list_of_old_polygons=[list_of_old_polygons;list_of_old_current_polygons];
        end
        % plot the resulting 2D RVE - single wall
%         for ii=1:1:size(list_of_old_elements,1)
%             line([list_of_old_elements(ii,2),list_of_old_elements(ii,4)],[list_of_old_elements(ii,3),list_of_old_elements(ii,5)],'Color','g');
%         end


        %% Create inner vertices
        % Preparing a list of vertices for mortar builder
        len=length(list_vertices_old_polygon);
        point_list_mortar_builder=list_vertices_old_polygon;
        
        % Delete the last entry of the list of vertices
        point_list_mortar_builder(len,:)=[];

        % Append the first entry of the vertex list to the last position.
        % And the last entry to the first position. So there is a list of
        % more than a circle of vertices
        point_list_mortar_builder_plus=[point_list_mortar_builder(length(point_list_mortar_builder),:);point_list_mortar_builder;point_list_mortar_builder(1,:)];
            
        % A List for the vertices of the new built inner polygon
        list_new_vertex=[];

        % Getting 3 connected points (equivalent 2 lines) to buid up
        % the current inner angle and the new vertex 'P_new' of the
        % inner polygon
        for h=2:1:length(point_list_mortar_builder_plus)-1

            % The 3 connected points are: 'P_previous', 'P_old_vertex', 'P_next'
            P_next = [point_list_mortar_builder_plus(h-1,1),point_list_mortar_builder_plus(h-1,2)];
            P_old_vertex = [point_list_mortar_builder_plus(h,1),point_list_mortar_builder_plus(h,2)];
            P_previous = [point_list_mortar_builder_plus(h+1,1),point_list_mortar_builder_plus(h+1,2)];
            % NOTE: Voronoi always builds clockwise polygons now

            % Calculate & plot the gradient of the first line
            gradient_vector_1 = [P_old_vertex(1,1)-P_previous(1,1);P_old_vertex(1,2)-P_previous(1,2)];
            gradient_1=gradient_vector_1(2,1)/gradient_vector_1(1,1);
%     			 line([P_old_vertex(1,1),P_previous(1,1)],[P_old_vertex(1,2),P_previous(1,2)],'Color','m')
%                hold on

            % Calculate & plot the gradient of the secound line
            gradient_vector_2 = [P_old_vertex(1,1)-P_next(1,1);P_old_vertex(1,2)-P_next(1,2)];
            gradient_2=gradient_vector_2(2,1)/gradient_vector_2(1,1);
%                 line([P_old_vertex(1,1),P_next(1,1)],[P_old_vertex(1,2),P_next(1,2)],'Color','g')
%                 hold on

            % Calculate the internal angle of the two lines
%                 internal_angle_2 = atan2(abs(det([P_next-P_old_vertex;P_previous-P_old_vertex])),dot(P_next-P_old_vertex,P_previous-P_old_vertex))
            internal_angle=acos((dot(gradient_vector_1,gradient_vector_2))...
                /(norm(gradient_vector_1)*norm(gradient_vector_2)));
            internal_angle_degree=internal_angle*180/pi;
            internal_angle=internal_angle/2;

            % Calculate the angle from origin (phi) to one current line of polygon 
            phi_old_rad = atan(gradient_1);
            phi_old_degree=phi_old_rad*180/pi;

            % Calculate the distance between old and new vertex
            hypotenuse=t_matrix/(2*sin(internal_angle));

            % Calculate the new angle from origin to new vertex
            phi_new_rad= phi_old_rad + internal_angle;
            phi_neu_degree=phi_new_rad*180/pi;

            % Calculate the displacment from old vertex to new vertex
            x_displacment = hypotenuse * cos(phi_new_rad);
            y_displacment = hypotenuse * sin(phi_new_rad);

            % Calculate new vertex 'P_new'
            P_new=[P_old_vertex(1,1) + x_displacment, P_old_vertex(1,2) + y_displacment];

            % Calculate whether P_new is inside the current polygon
            in = inpoly2(P_new,point_list_mortar_builder);

            % If it is not inside, mirror P_new
            if in==0
                P_new=[P_old_vertex(1,1) - x_displacment, P_old_vertex(1,2) - y_displacment];
            end
            
            % Write old and new vertex into a list, for some plot output
            list_old_to_new_vertex=[list_old_to_new_vertex;P_old_vertex,P_new];
            %distance_old_new = pdist(list_old_to_new_vertex);

            % Append the new vertex to the list of new vertices of the
            % current polygon
            list_new_vertex=[list_new_vertex;P_new(1,1),P_new(1,2)];
            
        end
        
        %disp('Create, cut and close inner polygons...')
        %% Call function for cutting and closing each inner cell at RVE border
        [list_of_current_elements,number_of_element,list_of_current_polygons,list_of_elements,]=...
            Cut_and_Close_a_Cell(list_new_vertex,list_of_current_elements,...
            number_of_element, area, L_RVE,list_of_elements,intersection_with_neighbour);
        
        % Handle the list_of_current_polygons to the global list_of_polygons
        if list_of_current_polygons ~= 0
            list_of_polygons=[list_of_polygons;list_of_current_polygons];
        end
    end
    
    %% Export data to input.txt for ABAQUS
    disp('Export data to ***.txt ...')
    % Collect all input parameters, which will be exported to python/ABAQUS
    input_file=[t_matrix;n_layer;L_RVE;t_cell;L_cell;t_matrix_layer;var_rnd_points/L_cell];
    file_name=strcat('_data\input.txt');
    file = fopen(file_name, 'w+');
    fprintf(file, '%i', input_file(1,1));
    fprintf(file, '\n');
    fprintf(file, '%i', input_file(2,1));
    fprintf(file, '\n');
    fprintf(file, '%i', input_file(3,1));
    fprintf(file, '\n');
    fprintf(file, '%i', input_file(4,1));
    fprintf(file, '\n');
    fprintf(file, '%i', input_file(5,1));
    fprintf(file, '\n');
    fprintf(file, '%i', input_file(6,1));
    fprintf(file, '\n');
    fprintf(file, '%i', input_file(7,1));
    fclose(file);
    
    % Export 'list_of_polygons' to list_of_polygons.txt for ABAQUS
    str_num_lay=num2str(num_lay);
    file_name=strcat('_data\list_of_polygons_', str_num_lay,'.txt');
    file = fopen(file_name, 'w+');

    for j=1:1:length(list_of_polygons)
        A=cell2mat(list_of_polygons(j,1));
        for i=1:1:size(A,2)
            if i~=size(A,2)
                fprintf(file, '%i\t', A(1,i));
            else
                fprintf(file, '%i', A(1,i));
            end    
        end
        if j~=length(list_of_polygons)
            fprintf(file, '\n');
        end 
    end
    fclose(file);

    % Export 'list_of_elements' to list_of_elements.txt for ABAQUS
    file_name=strcat('_data\list_of_elements_', str_num_lay,'.txt');
    file = fopen(file_name, 'w+');
    for j=1:1:length(list_of_elements)
        for i=1:1:size(list_of_elements,2)
            if i~=size(list_of_elements,2)
                fprintf(file, '%i\t', list_of_elements(j,i));
            else
                fprintf(file, '%i', list_of_elements(j,i));
            end    
        end
        if j~=length(list_of_elements)
            fprintf(file, '\n');
        end 
    end
    fclose(file);
    
    % Export 'list_of_old_polygons' to list_of_old_polygons.txt for ABAQUS
    str_num_lay=num2str(num_lay);
    file_name=strcat('_data\list_of_old_polygons_', str_num_lay,'.txt');
    file = fopen(file_name, 'w+');
    for j=1:1:length(list_of_old_polygons)
        A=list_of_old_polygons{j};
        for i=1:1:size(A,2)
            if i~=size(A,2)
                fprintf(file, '%i\t', A(1,i));
            else
                fprintf(file, '%i', A(1,i));
            end    
        end
        if j~=length(list_of_old_polygons)
            fprintf(file, '\n');
        end 
    end
    fclose(file);

    % Export 'list_of_old_elements' to list_of_old_elements.txt for ABAQUS

    file_name=strcat('_data\list_of_old_elements_', str_num_lay,'.txt');
    file = fopen(file_name, 'w+');
    for j=1:1:length(list_of_old_elements)
        for i=1:1:size(list_of_old_elements,2)
            if i~=size(list_of_old_elements,2)
                fprintf(file, '%i\t', list_of_old_elements(j,i));
            else
                fprintf(file, '%i', list_of_old_elements(j,i));
            end    
        end
        if j~=length(list_of_old_elements)
            fprintf(file, '\n');
        end 
    end
    fclose(file);

    %% Plot-Options the 2D-RVE

    %Size of points in some plots
    sz=100;
    %Line width for some plots
    grey=[0.7, 0.7, 0.7];
    blue=[49/255, 4/255, 180/255];
    green=[125/255, 205/255, 35/255];
    lw=3;
    lw_border=3;
    
    if graph~=0
        disp('Plot 2D RVE ...')
        
        % Set the scale of the plot axis
        xlim([scale_low_RVE scale_high_RVE])
        ylim([scale_low_RVE scale_high_RVE])
        axis equal 
        title('Periodic 2D RVE')
    end
    
        %% Plot the resulting 2D RVE
    for i=1:1:size(list_of_elements,1)
        line([list_of_elements(i,2),list_of_elements(i,4)],[list_of_elements(i,3),list_of_elements(i,5)],'LineWidth',lw,'Color',blue)
        hold on
    end

    if graph==1
        set(f1,'Visible','on')
        
        % Plot single-wall polygons in grey
        for i = 1:1:length(list_of_old_elements)
            line([list_of_old_elements(i,2), list_of_old_elements(i,4)],...
                [list_of_old_elements(i,3), list_of_old_elements(i,5)]...
                ,'LineWidth',lw,'Color',grey)
            hold on
        end
        
        % Plot line between old and new vertices in green
        for i = 1:1:length(list_old_to_new_vertex)
            line([list_old_to_new_vertex(i,1) list_old_to_new_vertex(i,3)],...
                [list_old_to_new_vertex(i,2) list_old_to_new_vertex(i,4)],'LineWidth',lw,'Color','k','LineStyle',':');
        end
        
        for i = 1:1:length(list_old_to_new_vertex)
        % Plot the new/inner vertices in blue
            scatter(list_old_to_new_vertex(i,3),list_old_to_new_vertex(i,4),sz,...
                'filled','MarkerEdgeColor',blue,'MarkerFaceColor',blue)
        end
        
        % Plot the old/single-wall vertices in grey    
        for i = 1:1:length(list_old_to_new_vertex)
            scatter(list_old_to_new_vertex(i,1),list_old_to_new_vertex(i,2),sz,...
                'filled','MarkerEdgeColor',grey,'MarkerFaceColor',grey)
        end
        
        % Plot all rnd_points
        for i = 1:1:size(periodic_point_list)
            scatter(periodic_point_list(i,1),periodic_point_list(i,2),sz,'filled','m')
        end
        
        % Plot the rnd_points inside the RVE
        for i = 1:1:size(list_of_rnd_points)
            scatter(list_of_rnd_points(i,1),list_of_rnd_points(i,2),sz,'filled','k')
        end
        % Plot the RVE+ area borders in black and ---
        area_2=[-L_RVE,-L_RVE; -L_RVE,2*L_RVE; 2*L_RVE,2*L_RVE; 2*L_RVE,-L_RVE; -L_RVE,-L_RVE];
        area_3=[0,-L_RVE;0,2*L_RVE;L_RVE,2*L_RVE;L_RVE,-L_RVE;0,-L_RVE];
        area_4=[-L_RVE,0;-L_RVE,L_RVE;2*L_RVE,L_RVE;2*L_RVE,0;-L_RVE,0];

        
        plot(area_2(:,1),area_2(:,2),'LineWidth',lw_border,'Color','k','LineStyle','--')
        plot(area_3(:,1),area_3(:,2),'LineWidth',lw_border,'Color','k','LineStyle','--')
        plot(area_4(:,1),area_4(:,2),'LineWidth',lw_border,'Color','k','LineStyle','--')
        plot(area(:,1),area(:,2),'LineWidth',lw_border,'Color','k')       
     
    end
    
    
    fprintf('*** End of Layer %d / %d *** \n',num_lay,n_layer)
end

%% End of program
disp('------------End of program------------')
