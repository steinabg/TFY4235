% Task 2.1e
% Calculate the sheet-to-sheet distance dxi_n for different dx and dt
% Input: s(x,t) vector
% Output: a vector = vector of distances between neighbour peaks in s
function a = task21e(s)

jmax = length(s);
[values, locations] = findpeaks(s); % Find local maxima in s
xspace = linspace(0,1,jmax);

a = zeros(length(values)-1,1);
for i = 1:length(values)-1
    a(i) = xspace(locations(i+1)) - xspace(locations(i));
end