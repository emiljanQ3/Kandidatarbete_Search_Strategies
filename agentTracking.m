function M = agentTracking(file)
    % Imports an xml-file
    % Converts the file to a matrix:
    % 1st dimenstion is the time length of the tracking
    % 2nd dimension is the number of parameters for each agent (t,x,y,z)
    % 3rd dimension is the number of agents (found by the tracking
    % programme)
    allTracks = importTrackMateTracks(file);
    length(allTracks);
    max_t = 1;
    for i = 1:length(allTracks)
        if (max(size(allTracks{i})) > max_t)
            max_t = max(size(allTracks{i}));
        end
    
    M = zeros(max_t,4,length(allTracks));
    lengths = zeros(length(allTracks));
    hold on    
    for i = 1:length(allTracks)
        m = cell2mat(allTracks(i));
        lengths(i) = size(m,1);
        M(1:size(m,1),:,i) = m;
        plot(m(:,2),m(:,3))
    end
    %lengths(18)
    %plot(M(1:lengths(18),2,18),M(1:lengths(18),3,18))
end

