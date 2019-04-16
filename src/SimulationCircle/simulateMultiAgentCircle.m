function [pos_a, rot_a, colision, totalTimeSteps] = simulateMultiAgentCircle(T,dT,D_r,D_p,v,W,numAgents,pos_a, rot_a, colision, R, r_c, interupt)

    totalTimeSteps = floor(T/dT);
    
    for T_i = 2:floor(T/dT)
        rot_a = mod(rot_a + dT * W' + sqrt(2 * D_r * dT) * randn(size(rot_a)), 2  * pi); %Update agent rotation for all agents

        targetPos = pos_a(:, :, T_i-1) + [cos(rot_a), sin(rot_a)] * dT * v + randn(numAgents, 2) * sqrt(2 * D_p * dT); %Calculate where a unhindered move would go.

        [pos_a(:, :, T_i), rot_a, colision(:,T_i)]= moveAllAgentsCircle(pos_a(:, :, T_i-1), targetPos,rot_a, R, r_c);    %Move agent and take obstacles into consideration.

        if(sum(colision(:,T_i)) && interupt)
            totalTimeSteps = T_i;
            return
        end
    end 
        
end