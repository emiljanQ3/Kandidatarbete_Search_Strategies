% Generates a labyrinth represented of size size(1)*size(j)
% At every index (i,j) in the labyrinth there is a a horizontal and a
% vertical line with probability p
% 1st dimension: x values of points in connected lines
% 2nd dimension: y values of points in connected lines
% 3rd dimension: 1 is horizontal lines, 2 is vertical lines
% 4th dimension: x position in labyrinth
% 5th dimension: y position in labyrinth
function labyrinth = generateLabyrinth(p,size)
clf;


labyrinth  = zeros(2,2,2,size(1)-1,size(2)-1);

for i = 1:size(1)-1
    for j = 1:size(2)-1
        walls = binornd(1,p,[1,2]);
        if (walls(1) == 1)
            labyrinth(:,:,1,i,j) = [i j; i j+1];
        else
            labyrinth(:,:,1,i,j) = [NaN NaN; NaN NaN];
        end
        if (walls(2) == 1)
            labyrinth(:,:,2,i,j) = [i j; i+1 j];
        else
            labyrinth(:,:,2,i,j) = [NaN NaN; NaN NaN];
        end
    end
end

% Plot labyrinth
% hold on
% for i = 1:size(1)-1
%     for j = 1:size(2)-1
%         for k = 1:2
%             plot(labyrinth(:,1,k,i,j),labyrinth(:,2,k,i,j),'k');
%         end
%     end
end
