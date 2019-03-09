function [w] = getChiralitySpiral(pos_a,dT, stepSizeThreshold,n)

    str = squeeze(pos_a(1,1,:)) ~= 0;
    index = strfind(str', [1 0]);
    pos = pos_a(:,:,1:index);

    pos = pos_a;

    for i = 1:size(pos, 3)-n
        i
        %find segment
        segment = pos(:,:,i:(i+n));
        %translate to origin
        segment = segment - segment(:,:,1);
        
        %Find rotation angle
        dX = pos(:,1,i+1)-pos(:,1,i);
        dY = pos(:,2,i+1)-pos(:,2,i);
        angle = -atan2(dY,dX);

        %Rotate segment
        for t = 1:size(segment,3)
            segment(:,:,t) = [cos(angle), -sin(angle); sin(angle), cos(angle)] * segment(:,:,t)';
        end
        
        dX2 = segment(:,1,2)-segment(:,1,1);
        dY2 = segment(:,2,2)-segment(:,2,1);
        angle2(i) = -atan2(dY2,dX2);
        
        %Debug below
        hold on
        plot(squeeze(segment(1,1,:)),squeeze(segment(1,2,:)))
        axis equal
    end
    w = 5
end