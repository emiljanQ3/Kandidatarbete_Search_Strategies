function [w] = getChiralitySpiral(pos,dT, stepSizeThreshold,n)

    for i = 2:size(pos, 2)
        %find segment
        segment = pos(:,i:i+n);
        %translate to origin
        segment = segment - segment(:,1);
        %Find rotation angle
        dX = pos(1,i+1)-pos(1,i-1);
        dY = pos(2,i+1)-pos(2,i-1);
        angle = tan(dY/dX);
        %Rotate segment
        segment = segment * [cos(angle), -sin(angle); sin(angle), cos(angle)];
        %Debug below
        plot(segment(1,:),segment(2,:));
    end

end