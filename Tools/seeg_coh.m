% This function takes a single matrix containing channels (rows) * time
% (columns) and calculates the channel-to-channel coherence between them,
% resulting in a single coherence based adjacency matrix

function coh = seeg_coh(M, ds)
if nargin < 2, ds = 1; end
coh = zeros(size(M,1));
count = 0;
textprogressbar('Calculating ');
for c = 1:size(M,1)
for cc = 1:size(M,1)
    count = count + 1;
    textprogressbar(100 * count / size(M,1)^2);
    m         = resample(M(c,:), 1, ds);
    mm        = resample(M(cc,:), 1, ds);
    coh(c,cc) = mean(mscohere(m, mm));
end
end
textprogressbar(' Done');