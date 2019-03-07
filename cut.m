function cutM = cut(file,agent)
% Reads an xml file into a matrix, for a specific agent
% Removes all NaNs in the beginning and end of the matrix
[M,lengths,times] = agentTracking(file);
ind = 1;
cutM = M(agent,:,times(agent,1):times(agent,2));
end 
