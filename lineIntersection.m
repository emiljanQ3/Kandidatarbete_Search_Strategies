function P = lineIntersection(line1, line2)
% line1 is the path of agent
% line2 is a wall

G = [5 6];


A = [line1(2,1)-line1(1,1) line2(1,1)-line2(2,1);
    line1(2,2)-line1(1,2) line2(1,2)-line2(2,2)];

B = [line2(1,1)-line1(1,1);
    line2(1,2)-line1(1,2)];

X = A\B;


% Find

if (X(1) > 1 || X(1) < 0 || X(2) > 1 || X(2) < 0)
    P = Inf*0;
    
    % TODO: Implement case when lines are "parallel" and intersect
    
    %     if (zeroDist)
    %
    %         %P = closest(line2)
    %         d1 = pdist([line1(1,:),line2(1,:)]);
    %         d2 = pdist([line1(1,:),line2(2,:)]);
    %         d3 = pdist([line1(2,:),line2(1,:)]);
    %         d4 = pdist([line1(2,:),line2(2,:)]);
    %
    %         if (d1 + d2 > d3 + d4)
    %             if (d1 < d2)
    %                 P = line2(1,:);
    %             else
    %                 P = line2(2,:);
    %             end
    %         else
    %             if (d3 < d4)
    %                 P = line2(1,:);
    %             else
    %                 P = line2(2,:);
    %             end
    %         end
    
else
    x = line1(1,1)+ X(1)*(line1(2,1)-line1(1,1));
    y = line1(1,2) + X(1)*(line1(2,2)-line1(1,2));
    
    P = [x,y];
end
end