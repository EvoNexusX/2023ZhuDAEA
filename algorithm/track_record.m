function [track] = track_record(groups)
    if ~isempty(groups)
        track = single([groups.xmean]');
    else
        track = single([]);
    end
end