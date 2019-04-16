function [pos_a,lengths,times] = agentTracking(file)
    % Imports an xml-file
    % Converts the file to a matrix:
    % 1st dimenstion is the number of agents (found by the tracking
    % programme)is the time length of the tracking
    % 2nd dimension is the x and y value
    % 3rd dimension is the number of time steps
    % Also returns the lengths

    allTracks = importTrackMateTracks(file);
    max_t = 1;
    for agent = 1:length(allTracks)
        if (size(allTracks{agent}, 1) > max_t)
            max_t = size(allTracks{agent},1);
        end
    end
    
    
    lengths = zeros(length(allTracks),1);
    pos_a = zeros(length(allTracks),2,max_t);
    times = zeros(size(pos_a,1),2);
    
    for agent = 1:length(allTracks)
        m = cell2mat(allTracks(agent));
        t_start = m(1,1)+1;
        t_end   = t_start+size(m,1)-1;
        
        lengths(agent) = size(m,1); % Actual number of points registered, not difference in time
        m_temp = m(:,2:3)';
        
        pos_a(agent,:,t_start:t_end) = m_temp;
        pos_a(agent,:,1:t_start-1) = NaN;
        pos_a(agent,:,t_start+size(m,1):end) = NaN;
        
        times(agent, 1) = t_start;
        times(agent, 2) = t_end;
    end
    
    
end

