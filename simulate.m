function [pos_a, rot_a, colision] = simulate(T,dT,D_r,D_p,v,w_i,numAgents,pos_a, rot_a,colision, obstacle, L, r_c,mapSize)

        for T_i = 2:floor(T/dT)
            rot_a = mod(rot_a + dT * w_i + sqrt(2 * D_r * dT) * randn(size(rot_a)), 2  * pi); %Update agent rotation for all agents
            targetPos = pos_a(:, :, T_i-1) + [cos(rot_a), sin(rot_a)] * dT * v + randn(numAgents, 2) * sqrt(2 * D_p * dT); %Calculate where a unhindered move would go.
            [pos_a(:, :, T_i), rot_a, colision(:,T_i)]= moveAllAgents(pos_a(:, :, T_i-1), targetPos,rot_a, obstacle, L, v*dT/10, r_c, mapSize);    %Move agent and take obstacles into consideration.        
        end 
        
        
end