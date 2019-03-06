function [pos_a, rot_a, colision, totalTime] = simulate(T,dT,D_r,D_p,v,w_i,numAgents,pos_a, rot_a,colision, obstacle, L, r_c,mapSize,edge)


        minPos = -floor(mapSize/2);
        totalTime = T;
        for T_i = 2:floor(T/dT)
            rot_a = mod(rot_a + dT * w_i + sqrt(2 * D_r * dT) * randn(size(rot_a)), 2  * pi); %Update agent rotation for all agents
            targetPos = pos_a(:, :, T_i-1) + [cos(rot_a), sin(rot_a)] * dT * v + randn(numAgents, 2) * sqrt(2 * D_p * dT); %Calculate where a unhindered move would go.
            [pos_a(:, :, T_i), rot_a, colision(:,T_i)]= moveAllAgents(pos_a(:, :, T_i-1), targetPos,rot_a, obstacle, L, v*dT/10, r_c, mapSize);    %Move agent and take obstacles into consideration.        

            if(edge)
                cellIndex = floor((pos_a(:,:,T_i) - minPos)/L) + 1;
                if any(cellIndex(1) > mapSize(1) || cellIndex(2) > mapSize(2) || cellIndex(1) < 1 || cellIndex(2) < 1)
                   
                   pos_a = pos_a(:,:,1:T_i);
                   totalTime = (T_i-1)*dT;
                   break
                end
            end
        end 

        
end