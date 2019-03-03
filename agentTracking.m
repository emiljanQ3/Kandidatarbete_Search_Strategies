function [M,lengths,times] = agentTracking(file)
    % Imports an xml-file
    % Converts the file to a matrix:
    % 1st dimenstion is the number of agents (found by the tracking
    % programme)is the time length of the tracking
    % 2nd dimension is the x and y value
    % 3rd dimension is the number of time steps
    % Also returns the lengths

    allTracks = importTrackMateTracks(file);
    length(allTracks);
    max_t = 1;
    for i = 1:length(allTracks)
        if (max(size(allTracks{i})) > max_t)
            max_t = max(size(allTracks{i}));
        end
    max_t;
    M = zeros(length(allTracks),2,max_t);
    lengths = zeros(length(allTracks),1);
    hold on    
    for i = 1:length(allTracks)
        i;
        m = cell2mat(allTracks(i));
        t_start = m(1,1)+1;
        lengths(i) = size(m,1); % Actual number of points registered, not difference in time
        m_temp = m(:,2:3)';
        M(i,:,t_start:t_start+size(m,1)-1) = m_temp;
        M(i,:,1:t_start-1) = NaN;
        M(i,:,t_start+size(m,1):end) = NaN;
        t_start+size(m,1);
        max_t;
        %plot(m(:,2),m(:,3))
    end
    
    times = zeros(size(M,1),2);
    for agent = 1: size(M,1)
        agent;
        t = 1;
        while (isnan(M(agent,:,t)))
            t=t+1;
        end
        times(agent,1) = t; % Start time
        end_time = t+lengths(agent)-1;
        end_time;
        times(agent,2) = t+lengths(agent)-1; % End time
    end
%     %lengths(18)
%     X_vals = M(3,1,:);
%     Y_vals = M(3,2,:);
%     X_vals = X_vals(:,:)';
%     Y_vals = Y_vals(:,:)';
%     plot(X_vals,Y_vals)
lengths;
end

