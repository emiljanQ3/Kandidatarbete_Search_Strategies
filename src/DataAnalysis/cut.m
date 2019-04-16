function [cutPos_a,length, times] = cut(file,agent)
% Reads an xml file into a matrix, for a specific agent
% Removes all NaNs in the beginning and end of the matrix
[pos_a,lengths,times] = agentTracking(file);
cutPos_a = pos_a(agent,:,times(agent,1):times(agent,2));
length = lengths(agent);
times = times(agent,:);
end 
