function loc = artifact_detect(window, threshold)
% ARTIFACT DETECT takes in a window and returns the location of an
% artificial pacing artifact, if present. If none is present it returns NaN.
% A peak is an artifact if it is greater than threshold.
% loc: the location if present. If not present, loc is NaN.
% window: a window of the signal to be processed.
% threshold: the threshold which is used for artifact detection. This is
% an optional argument

%     if nargin == 1
%         threshold = 2000; % TODO come up with good default value
%     end
    
    abs_window = abs(window);
    [peak_mag, peak_loc] = max(abs_window);
    
    if peak_mag > threshold
        loc = peak_loc;
    else
        loc = NaN;
    end
    
end


