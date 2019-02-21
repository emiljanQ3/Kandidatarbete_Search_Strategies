function [M,lengths] = agentTracking(file)
    % Imports an xml-file
    % Converts the file to a matrix:
    % 1st dimenstion is the number of agents (found by the tracking
    % programme)is the time length of the tracking
    % 2nd dimension is the time length of the trackingis the number of parameters for each agent (t,x,y,z)
    % 3rd dimension is the number of parameters for each agent (t,x,y,z)

    allTracks = importTrackMateTracks(file);
    length(allTracks);
    max_t = 1;
    for i = 1:length(allTracks)
        if (max(size(allTracks{i})) > max_t)
            max_t = max(size(allTracks{i}));
        end
    
    M = zeros(length(allTracks),2,max_t);
    lengths = zeros(length(allTracks),1);
    hold on    
    for i = 1:length(allTracks)
        m = cell2mat(allTracks(i));
        lengths(i) = size(m,1);
        m_temp = m(:,2:3)';
        M(i,:,1:size(m,1)) = m_temp;
        plot(m(:,2),m(:,3))
    end
    %lengths(18)
    %plot(M(18,1:lengths(18),2),M(18,1:lengths(18),3))
end

