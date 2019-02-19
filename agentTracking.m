function M = agentTracking(file)

    allTracks = importTrackMateTracks(file);
    length(allTracks);
    M = zeros(893,4,length(allTracks));
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

