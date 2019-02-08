function newPos = moveAgent(pos,rot,dist,obstacle,L)
         targetPos = [cos(rot),sin(rot)]*dist;
         newPos = moveAgentPrime(pos, targetPos, obstacle, L, dist/100);
         
         
end