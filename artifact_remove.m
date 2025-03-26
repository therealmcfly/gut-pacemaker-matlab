function output = artifact_remove(window, x1, x2)
%ARTIFACT REMOVE takes in a window of a signal containing a pacing artifact
%starting at location x and replaces it with a hyperbolic cosine function.
%   WINDOW: a window of the signal
%   X1: the index of the sample at which the artifact starts
%   X2: the index of the sample at which the artifact ends
%   OUTPUT: the same window of the signal, now with the artifact removed

% compute start and end slopes
dy1 = window(x1) - window(x1-1);
dy2 = window(x2) - window(x2-1);

% calculate spline with sample rate at 32 Hz
x = [x1 x2];
y = [window(x1) window(x2)];
cs = spline(x,[dy1 y dy2]); % cubic spline
%xx = x1:1/32:x2;
xx = x1:x2;
replace = ppval(cs,xx);

% output result
output = zeros(1,length(window));
output(1:x1-1) = window(1:x1-1);
output(x1:x2) = replace;
output(x2+1:end) = window(x2+1:end);

end
