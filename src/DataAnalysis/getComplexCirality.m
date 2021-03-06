function [w_average, v] = getComplexCirality(pos_a, dT, stepSizeThreshold)
   
   
    w_a = zeros(1,length(pos_a(:,1,1)));
    v_a = zeros(1,length(pos_a(:,1,1)));
    n_a = zeros(1,length(pos_a(:,1,1)));
    
    for agent = 1:length(pos_a(:,1,1))
        
        str = squeeze(pos_a(agent,1,:)) ~= 0;
        index = strfind(str', [1 0]);
        temp_pos = pos_a(agent,:,1:index);
        
        [w_a(agent), v_a(agent),wigth(agent)] = getCirality(temp_pos, dT, stepSizeThreshold);
        n_a(agent) = length(temp_pos(1,1,:)); 

    end

    w_average = sum(w_a.*wigth)/sum(wigth);
    v = mean(v_a);
end


