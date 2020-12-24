function [tracks] = importTrackImaris(file_path)
%% file_path: the path to the csv file from Imaris statistic position output of tracking module 
    T = readtable(file_path);
    track_IDs = unique(T.TrackID);
    tracks = {};
    for t = 1:length(track_IDs)
        tracks(t) = { T{T.TrackID==track_IDs(t),{'Time', 'PositionX', 'PositionY', 'PositionZ'}} };
    end
    
end