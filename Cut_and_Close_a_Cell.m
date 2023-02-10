function [list_of_current_elements,number_of_element, list_of_current_polygons,...
    list_of_elements]=Cut_and_Close_a_Cell(list_vertexes_one_cell,list_of_current_elements,...
    number_of_element, area, L_RVE,list_of_elements,intersection_with_neighbour);

%% Check location of each line of current inner_polygon

% Add the last vertex at first position of the list
list_new_vertex_added=[list_vertexes_one_cell(length(list_vertexes_one_cell),:);list_vertexes_one_cell];
list_of_polygons={};
% Pick two nearby points (= one element) out of the inner_polygon
for m=1:1:length(list_new_vertex_added)-1
    x1=list_new_vertex_added(m,1);
    y1=list_new_vertex_added(m,2);
    x2=list_new_vertex_added(m+1,1);
    y2=list_new_vertex_added(m+1,2);

    % Current inner_vertex_line
    L1=[x1 x2;y1 y2];

    %% Check wheather a line is completly inside the RVE
    if x1>=0 && y1>=0 && x1<=L_RVE && y1<=L_RVE && x2>=0 && y2>=0 && x2<=L_RVE && y2<=L_RVE
        x=[x1,y1;x2,y2];
        % Plot all lines, which are completly inside the RVE in
        % red
        %line(x(:,1),x(:,2),'LineWidth',lw_rve,'Color','r')
        %hold on

        number_of_element=number_of_element+1;
        new_element=[number_of_element,x1,y1,x2,y2];
        list_of_current_elements=[list_of_current_elements;new_element];
    end

    %% Check intersection(s) with RVE border
    double_intersection=[];
    z=0;
    kk=0;
    % Check intersection(s) for every (value=4) border_line seperated
    for k=1:1:4
        x3=area(k,1);
        y3=area(k,2);
        x4=area(k+1,1);
        y4=area(k+1,2);

        % Current border_line
        L2=[x3 x4;y3 y4];

        % Check each intersection between current inner_polygon_line and current border_line
        P = InterX(L1,L2);
        TF = isempty(P);

        % If there is a intersection, then...
        if TF==0
            % Remember the No. of border line, which gets
            % intersected.
            intersection_border=k;
            if k==1
                kk=1;
            end
            % The intersection function 'interX' shows
            % sometimes a mistake in the calulation of 10^-15. This is my
            % solution:
            Differenz_1=abs(P(1,1)-L_RVE);
            Differenz_2=abs(P(2,1)-L_RVE);
            if Differenz_1 < 0.00000001
                P(1,1)=round(P(1,1)*100000000)/100000000;
            elseif Differenz_2 < 0.00000001
                P(2,1)=round(P(2,1)*100000000)/100000000;
            end    
            % Counter for intersections per line & save the
            % intersection
            z=z+1;
            double_intersection=[double_intersection;P(1,1),P(2,1)];

            % if there are two intersections with a borders,
            % then...
            if z==2
                xx2=double_intersection(1,1);
                yy2=double_intersection(1,2);
                xx1=double_intersection(2,1);
                yy1=double_intersection(2,2);

%                 % Plot the two-intersection-lines
%                 line(double_intersection(:,1),double_intersection(:,2),'LineWidth',lw_rve,'Color','m');
%                 hold on

                % Add the ploted line to output-list: 'list_of_current_elements'
                number_of_element=number_of_element+1;

                % Check whether both intersections points are in the right direction
                if k==4 && kk==1
                    new_element=[number_of_element,xx2,yy2,xx1,yy1];
                else
                    new_element=[number_of_element,xx1,yy1,xx2,yy2];
                end
                list_of_current_elements=[list_of_current_elements;new_element];

                % Plot the two intersection-points with RVE border, in magenta
%                 sz=15;
%                 scatter(double_intersection(:,1),double_intersection(:,2),sz,'filled','g')
            end
        end

        % If there is only one intersection with RVE border
        if k==4 && z ==1
        x_intersect=double_intersection(1,1);
        y_intersect=double_intersection(1,2);
            % Prove, if the original element starts/ends at RVE
            % Border
            if x1~=0 && x1~=L_RVE && y1~=0 && y1~=L_RVE && x2~=0 && x2~=L_RVE && y2~=0 && y2~=L_RVE
                    % Plot intersection-point at RVE border, in blue
                    %sz=15;
                    %scatter(x_intersect,y_intersect,sz,'filled','b')

                    % Determine the direction of intersection-line
                    % If the intersection-line intersect first border,
                    % then...
                if intersection_border==1
                    if x1>=0
%                         line([x1,x_intersect],[y1,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx1=x1;
                        yy1=y1;
                        xx2=x_intersect;
                        yy2=y_intersect;

                    elseif x2>=0
%                         line([x2,x_intersect],[y2,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx2=x2;
                        yy2=y2;
                        xx1=x_intersect;
                        yy1=y_intersect;
                    end
                % If the intersection-line intersect secound border,
                % then...
                elseif intersection_border==2
                    if y1<=L_RVE
%                         line([x1,x_intersect],[y1,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx1=x1;
                        yy1=y1;
                        xx2=x_intersect;
                        yy2=y_intersect;
                    elseif y2<=L_RVE
%                         line([x2,x_intersect],[y2,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx2=x2;
                        yy2=y2;
                        xx1=x_intersect;
                        yy1=y_intersect;
                    end
                % If the intersection-line intersect third border,
                % then...
                elseif intersection_border==3
                    if x1<=L_RVE
%                         line([x1,x_intersect],[y1,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx1=x1;
                        yy1=y1;
                        xx2=x_intersect;
                        yy2=y_intersect;
                    elseif x2<=L_RVE
%                         line([x2,x_intersect],[y2,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx2=x2;
                        yy2=y2;
                        xx1=x_intersect;
                        yy1=y_intersect;
                    end
                % If the intersection-line intersect fourth border,
                % then...
                elseif intersection_border==4
                    if y1>=0
%                         line([x1,x_intersect],[y1,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx1=x1;
                        yy1=y1;
                        xx2=x_intersect;
                        yy2=y_intersect;                            
                    elseif y2>=0
%                         line([x2,x_intersect],[y2,y_intersect],'LineWidth',lw_rve,'Color','r')
%                         hold on
                        xx2=x2;
                        yy2=y2;
                        xx1=x_intersect;
                        yy1=y_intersect;
                    end
                end

                % Add the ploted line to output-list: 'list_of_current_elements'
                number_of_element=number_of_element+1;
                new_element=[number_of_element,xx1,yy1,xx2,yy2];
                list_of_current_elements=[list_of_current_elements;new_element];
            end
        end  
    end

    %% Close the cells which are devided by  RVE border
    % If the current inner polygon has more than 2 elements inside RVE-area, then...
    if size(list_of_current_elements,1) > 1
        last_element=size(list_of_current_elements,1);
        next_to_last_element=size(list_of_current_elements,1)-1;

        % Prove the internal structur data of the saved points.
        % If this is true, then the intersections hit two
        % different borders
        if list_of_current_elements(next_to_last_element,4) ~= list_of_current_elements(last_element,2) && list_of_current_elements(next_to_last_element,5) ~= list_of_current_elements(last_element,3) && list_of_current_elements(next_to_last_element,4) ~= list_of_current_elements(last_element,4) && list_of_current_elements(next_to_last_element,5) ~= list_of_current_elements(last_element,5)
            %disp('There are two points not on the same border. But there is no Problem.')

            % Check corner intersection and connect those
            xx1=list_of_current_elements(next_to_last_element,4);
            yy1=list_of_current_elements(next_to_last_element,5);
            xx3=list_of_current_elements(last_element,2);
            yy3=list_of_current_elements(last_element,3);

            % Prove which edge of the RVE-area, will be
            % connected to both intersection points
            if yy1 == 0 || yy1 == L_RVE  
                yy2=yy1;
            elseif xx1 == 0 || xx1 == L_RVE
                xx2=xx1;
            end
            if  yy3 == 0 || yy3 == L_RVE
                yy2=yy3;
            elseif xx3 == 0 || xx3 == L_RVE
                xx2=xx3;
            end  

            % Plot the connection to/from the corner, in cyan
%             line([xx1,xx2],[yy1,yy2],'LineWidth',lw_rve,'Color','c')
            number_of_element=number_of_element+1;
            new_element_1=[number_of_element,xx1,yy1,xx2,yy2];
            
%             line([xx2,xx3],[yy2,yy3],'LineWidth',lw_rve,'Color','c')
            number_of_element=number_of_element+1;
            new_element_2=[number_of_element,xx2,yy2,xx3,yy3];

            % Add the ploted lines to output-list: 'list_of_current_elements'
            list_of_current_elements=[list_of_current_elements;new_element_1;new_element_2];
            list_of_current_elements([size(list_of_current_elements,1)-2,size(list_of_current_elements,1)],[2,3,4,5]) = list_of_current_elements([size(list_of_current_elements,1),size(list_of_current_elements,1)-2],[2,3,4,5]);
            list_of_current_elements([size(list_of_current_elements,1)-1,size(list_of_current_elements,1)-2],[2,3,4,5]) = list_of_current_elements([size(list_of_current_elements,1)-2,size(list_of_current_elements,1)-1],[2,3,4,5]);

            % Prove the internal structur data of the saved points.
            % If this is true, then the intersections hit the same
            % border line
        elseif list_of_current_elements(next_to_last_element,4) == list_of_current_elements(last_element,2) || list_of_current_elements(next_to_last_element,5) == list_of_current_elements(last_element,3)
            xx1=list_of_current_elements(next_to_last_element,4);
            yy1=list_of_current_elements(next_to_last_element,5);
            xx2=list_of_current_elements(last_element,2);
            yy2=list_of_current_elements(last_element,3);

            % Plot the line, in green
            %line([xx1,xx2],[yy1,yy2],'LineWidth','5','Color','g')

            % Add the ploted line to output-list: 'list_of_current_elements'
            number_of_element=number_of_element+1;
            new_element=[number_of_element,xx1,yy1,xx2,yy2];
            list_of_current_elements=[list_of_current_elements;new_element];
            list_of_current_elements([size(list_of_current_elements,1)-1,size(list_of_current_elements,1)],[2,3,4,5]) = list_of_current_elements([size(list_of_current_elements,1),size(list_of_current_elements,1)-1],[2,3,4,5]);                       
        end
    end
end

%% Add the geometric informations to the goldbal output lists

% Check if the current inner polygone has any elements inside or
% intersecting RVE-area. True, if there are elements.
if isempty(list_of_current_elements)==0
    last_element=size(list_of_current_elements,1);

    % Check whether polygon is already closed (last vertex is equal
    % fist vertex. True, if it is not closed.
    if list_of_current_elements(1,2) ~= list_of_current_elements(last_element,4) || list_of_current_elements(1,3) ~= list_of_current_elements(last_element,5)

        % Check whether the not-closed polygon is intersected on only
        % one border line
        if list_of_current_elements(1,2)== list_of_current_elements(last_element,4) || list_of_current_elements(1,3)== list_of_current_elements(last_element,5)
            xx2=list_of_current_elements(1,2);
            yy2=list_of_current_elements(1,3);
            xx1=list_of_current_elements(last_element,4);
            yy1=list_of_current_elements(last_element,5);

%             % Plot the line, in yellow
%             line([xx1,xx2],[yy1,yy2],'Color','y')

            % Add the ploted line to output-list: 'list_of_current_elements'
            number_of_element=number_of_element+1;
            new_element=[number_of_element,xx1,yy1,xx2,yy2];
            list_of_current_elements=[list_of_current_elements;new_element];

            % Check and correct the order of list_of_current_elements
            % switch last and beforelast element
            if list_of_current_elements(1,4)~=list_of_current_elements(2,2) && list_of_current_elements(1,5)~=list_of_current_elements(2,3)
                fprintf('LOOK HERE, IN THIS CASE THE PART IS NECESSARY')
                list_of_current_elements([size(list_of_current_elements,1)-1,size(list_of_current_elements,1)],[2,3,4,5]) = list_of_current_elements([size(list_of_current_elements,1),size(list_of_current_elements,1)-1],[2,3,4,5]);
            end

        % Not-closed polygon is intersected by two border lines
        else                    
            % Check in which corner the intersections are
            xx3=list_of_current_elements(1,2);
            yy3=list_of_current_elements(1,3);
            xx1=list_of_current_elements(last_element,4);
            yy1=list_of_current_elements(last_element,5);

            if list_of_current_elements(1,3) == 0 || list_of_current_elements(1,3) == L_RVE  
                yy2=yy3;
            elseif list_of_current_elements(1,2) == 0 || list_of_current_elements(1,2) == L_RVE
                xx2=xx3;
            end
            if  list_of_current_elements(last_element,5) == 0 || list_of_current_elements(last_element,5) == L_RVE
                yy2=yy1;
            elseif list_of_current_elements(last_element,4) == 0 || list_of_current_elements(last_element,4) == L_RVE
                xx2=xx1;
            end  

%             % Plot the connection to/from the corner, in green
%             line([xx1,xx2],[yy1,yy2],'LineWidth',lw_rve,'Color','b')
            number_of_element=number_of_element+1;
            new_element_1=[number_of_element,xx1,yy1,xx2,yy2];

%             line([xx2,xx3],[yy2,yy3],'LineWidth',lw_rve,'Color','b')
            number_of_element=number_of_element+1;
            new_element_2=[number_of_element,xx2,yy2,xx3,yy3];

            % Add the ploted lines to output-list: 'list_of_current_elements'
            list_of_current_elements=[list_of_current_elements;new_element_1;new_element_2];


        end
    end

    % zero all elements, which are builded by the same points
    for iii=1:1:length(list_of_current_elements(:,1))    
        if list_of_current_elements(iii,2)==list_of_current_elements(iii,4) && list_of_current_elements(iii,3)==list_of_current_elements(iii,5)
            list_of_current_elements(iii,:)=0;
        end
    end
    % delete all rows with only 'zero'-elements
    list_of_current_elements( all(~list_of_current_elements,2), : ) = [];

    %% In some cases it is possible that the inner polygon lines are intersecting each other. In this case the intersection point will be the new vertex of the polygon:
    % Check, if there is an intersection of two elements, which
    % have the same neighbour element
    if length(list_of_current_elements(:,1))>3 && intersection_with_neighbour==0
        iii=0;
        while iii<length(list_of_current_elements(:,1))
            iii=iii+1;
            E1=[list_of_current_elements(iii,2) list_of_current_elements(iii,4);list_of_current_elements(iii,3) list_of_current_elements(iii,5)];
            % next to last round of while-loop
            if iii==length(list_of_current_elements(:,1))-1
                E2=[list_of_current_elements(iii+1,2) list_of_current_elements(iii+1,4);list_of_current_elements(iii+1,3) list_of_current_elements(iii+1,5)];
                E3=[list_of_current_elements(1,2) list_of_current_elements(1,4);list_of_current_elements(1,3) list_of_current_elements(1,5)];
            % the last round of while-loop
            elseif iii==length(list_of_current_elements(:,1))
                E2=[list_of_current_elements(1,2) list_of_current_elements(1,4);list_of_current_elements(1,3) list_of_current_elements(1,5)];
                E3=[list_of_current_elements(2,2) list_of_current_elements(2,4);list_of_current_elements(2,3) list_of_current_elements(2,5)];
            else
                E2=[list_of_current_elements(iii+1,2) list_of_current_elements(iii+1,4);list_of_current_elements(iii+1,3) list_of_current_elements(iii+1,5)];
                E3=[list_of_current_elements(iii+2,2) list_of_current_elements(iii+2,4);list_of_current_elements(iii+2,3) list_of_current_elements(iii+2,5)];
            end
            P = InterX(E1,E3);
            TF = isempty(P);

            % If there is a intersection between E1 and E3, then...
            if TF==0
                list_of_current_elements(iii,4)=P(1,1);
                list_of_current_elements(iii,5)=P(2,1);                        
                if iii==length(list_of_current_elements(:,1))-1
                    list_of_current_elements(1,2)=P(1,1);
                    list_of_current_elements(1,3)=P(2,1);
                    list_of_current_elements(iii+1,:)=0;
                elseif iii==length(list_of_current_elements(:,1))
                    list_of_current_elements(2,2)=P(1,1);
                    list_of_current_elements(2,3)=P(2,1);
                    list_of_current_elements(1,:)=0;
                else
                    list_of_current_elements(iii+2,2)=P(1,1);
                    list_of_current_elements(iii+2,3)=P(2,1);
                    list_of_current_elements(iii+1,:)=0;
                end
                if length(list_of_current_elements(:,1))>3
                    iii=0;
                else
                    iii=length(list_of_current_elements(:,1));
                end
            end

            % delete all rows with only 'zero'-elements
            list_of_current_elements( all(~list_of_current_elements,2), : ) = [];
        end
    end

    % Add all elements of the current polygon to a layer-global
    % list of elements
    list_of_elements=[list_of_elements;list_of_current_elements];

    % Wirte the position of elements for the current inner polygon
    % into one row of 'list_of_polygons'
    list_of_current_polygons=[];
    for ii=1:1:size(list_of_current_elements,1)
        list_of_current_polygons=[list_of_current_polygons,list_of_current_elements(ii,1)];
    end    
else
    % gives this parameter just a value for returning
    list_of_current_polygons=0;
end