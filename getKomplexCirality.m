function [cir, v] = getKomplexCirality(pos_a,dT)

    
    for agent = 1:length(pos_a(:,1,1))

        str = squeeze(pos_a(agent,1,:)) ~= 0;
        index = strfind(str', [1 0]);
        temp_pos = pos_a(agent,:,1:index);
        
        n = 0;
        t = 1;

        u = temp_pos(1,:,2) - temp_pos(1,:,1); % om första är 0 eller liten? måste fixas stort fel
        
        v(1) = norm(u)/dT;
        orientation(1,agent) = atan2(u(2),u(1));
        T(1) = 0;
        if(norm(u)<1)
            orientation(1,agent) = 0
        end

        for i = 2:length(temp_pos)-1

            u = temp_pos(1,:,i+1) - temp_pos(1,:,i);
            v(i) = norm(u)/dT;

            if (norm(u)>1)
                t = t+1
                T(t) = dT*(i-1);
                angel =  atan2(u(2),u(1))
                if((orientation(t-1,agent)-(angel+2*pi*n)) > 5)
                    n = n+1;

                elseif((orientation(t-1,agent) - (angel+2*pi*n)) < -5)
                    n = n-1;
                end
                orientation(t,agent) = angel + 2*n*pi;
            end

        end
    % Orientation is calculated now to get the cirality        
        w = 0;
        for j = 1:length(orientation(:,agent))-1
            if(orientation(j+1,agent)~= 0 && orientation(j,agent)~= 0)
                j
                cir2(j) = (orientation(j+1,agent)-orientation(j,agent))/(T(j+1)-T(j)); 
                w = w +1;
            end
        end

        cir(1,agent) = mean(cir2);
        cir(2,agent) = w; 

    end

    cir = sum(cir(1,:).*cir(2,:))/sum(cir(2,:))
    v = mean(v);
end


