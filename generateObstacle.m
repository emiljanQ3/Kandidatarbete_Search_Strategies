%R and r represent legth of obstacle elements and are expressed as
%fractions of cell side length.

%The return value obstacle is a list of pairs of points, where eeach pair
%represents a line segment. For example the lines L_1 and L_2 could be
%stored as:

% obstacle(:,:,1) = [x_1, y_1; x_2, y_2]

function obstacle = generateObstacle(string, R, r )
    
    center = [1/2, 1/2];

    if string == "c" %Cross
        obstacle(:,:,1) = [-R/2, 0; R/2, 0] + center;
        obstacle(:,:,2) = [0, -R/2; 0, R/2] + center;
       return 
    end
    
    if string == "hc" %Hooked cross
        %Cross
        obstacle(:,:,1) = [-R/2, 0; R/2, 0] + center;
        obstacle(:,:,2) = [0, -R/2; 0, R/2] + center;
        
        %Hooks
        obstacle(:,:,3) = [-R/2, -r/2; -R/2, r/2] + center;
        obstacle(:,:,4) = [R/2, -r/2; R/2, r/2] + center;
        obstacle(:,:,5) = [-r/2, R/2; r/2, R/2] + center;
        obstacle(:,:,6) = [-r/2, -R/2; r/2, -R/2] + center;
        return
    end
    
    if string == "hm" %Homogenous cell
        obstacle = null;
        return
    end 
    
end