clear all
D         = seeg_housekeeping;
Fbase     = D.Fbase;
Fscripts  = D.Fscripts;
Fdata     = D.Fdata;
Fanalysis = D.Fanalysis;
fs        = filesep;

files       = cellstr(spm_select('FPList', [Fanalysis fs 'Win_Coh'], '.mat$'));

% Load data file to extract header labels
%--------------------------------------------------------------------------
hdr     = ft_read_header([Fdata fs 'Awake.edf']);
chanlab = hdr.label(2:39);

% Load individual files into concatenated matrix
%==========================================================================
allcoh = [];
for f = 1:length(files)
    
load(files{f});
wcoh    = B.wcoh;   clear B;
wl(f)   = size(wcoh,1);
dyncoh  = [];

for c = 1:size(wcoh,1)
    thiscoh      = squeeze(wcoh(c,:,:));
    triid        = find(tril(ones(size(wcoh,2)),-1));     % lower triangle indices
    dyncoh(:,c)  = thiscoh(triid);
end

allcoh = [allcoh, dyncoh];
end

%% Calculate errors for increasing numbers of subnetworks
%--------------------------------------------------------------------------
for k = 1:10
    [W H i t r] = nmfnnls(allcoh, k);
    Err(k)      = sum(sum(sqrt((allcoh - W*H).^2))) / sum(sum(allcoh));
    R(k)        = r;
end

% Plot results
%--------------------------------------------------------------------------
figure(1)
    plot(Err), hold on
    plot([0 length(Err)], [0.1 0.1]);
    scatter(1:length(Err), Err, 'k.');
    % Labels
    title('Residual error');
    ylabel('Error');
    xlabel('Number of Networks');    
    % Settings
    ylim([0 Inf]);

    min_k = find(Err < 0.1);    min_k = min_k(1);

%% Matrix decomposition with optimal number of subnetworks
%==========================================================================
if isempty(min_k), min_k = 3; end
[W H]   = nmfnnls(allcoh, min_k);

figure(2)
clear A G B
for w = 1:size(W,2)
    A{w}    = seeg_untril(W(:,w));
    B{w}    = A{w} > mean(A{w}(:));
    G{w}    = graph(B{w}, chanlab);
    
    subplot(2,min_k,w)
    plot(G{w}, 'Layout', 'circle')
    axis square
    subplot(2,min_k,w+min_k)
    imagesc(A{w});
    axis square
end

%% Plot coherence dynamics and subgraph expression plots
%--------------------------------------------------------------------------
figure(3)
for i = 1:min_k
    lbl{i} = num2str(i);
end

subplot(3,1,1)
imagesc(allcoh), hold on
for w = 1:length(wl)
    x = sum(wl(1:w));
    plot([x x], [0 1000], 'color', 'w', 'linewidth', 3);
end

subplot(3,1,2)
plot(H'); hold on
legend(lbl);
for w = 1:length(wl)
    x = sum(wl(1:w));
    plot([x x], [0 1000], 'color', 'k', 'linewidth', 3);
end
ylim([0 max(max(H))]);
xlim([1 size(H,2)]);

subplot(3,1,3)
change      = abs(diff(H'));
abschange   = smooth(sum(change,2));
[val loc]   = findpeaks(abschange);
peakchg     = loc(val > 1.5);

plot(abschange); hold on
plot([1 length(abschange)], [1.5 1.5], 'k');
xlim([0 Inf])
title('Changes in network expression');

%%
pi = [2 7 8];
subplot(3,1,2)
hold on
for p = pi
    plot([peakchg(p)+10 peakchg(p)+10],[0 max(H(:))],'color', [.5 .5 .5]);
    plot([peakchg(p)-10 peakchg(p)-10],[0 max(H(:))],'color', [.5 .5 .5]);
end

