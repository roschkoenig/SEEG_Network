clear all
D         = seeg_housekeeping;
Fbase     = D.Fbase;
Fscripts  = D.Fscripts;
Fdata     = D.Fdata;
Fanalysis = D.Fanalysis;
fs        = filesep;
addpath([Fscripts fs 'nmfv1_4']);

% Load coherence file
%==========================================================================
files   = cellstr(spm_select('List', Fanalysis, '^*win.mat'));
f       = 1;
plotem  = 1;

for f = 1:length(files)
    
load([Fanalysis fs files{f}]);
wcoh    = B.wcoh;
clear B;

for c = 1:size(wcoh,1)
    thiscoh      = squeeze(wcoh(c,:,:));
    triid        = find(tril(ones(size(wcoh,2)),-1));     % lower triangle indices
    dyncoh(:,c)  = thiscoh(triid);
end

% Apply matrix decomposition with sparsity constraints
%==========================================================================
[W H] = nmfnnls(dyncoh, 3)

if plotem
    figure(f)
    subplot(3,3,1:6)
        plot(H');
        legend({'First', 'Second', 'Third'});

    for p = 1:3
        subplot(3,3,6+p)
        imagesc(seeg_untril(W(:,p)));
        colorbar
        axis square
    end
end
end

%% Load individual files into concatenated matrix
%==========================================================================
allcoh = [];
for f = 1:length(files)
    
load([Fanalysis fs files{f}]);
wcoh    = B.wcoh;
clear B;
wl(f)   = size(wcoh,1)
dyncoh  = [];

for c = 1:size(wcoh,1)
    thiscoh      = squeeze(wcoh(c,:,:));
    triid        = find(tril(ones(size(wcoh,2)),-1));     % lower triangle indices
    dyncoh(:,c)  = thiscoh(triid);
end

allcoh = [allcoh, dyncoh];
end

figure 
lbl = {'1', '2', '3', '4'};
subplot(3,4,[1:4])
imagesc(allcoh), hold on
for w = 1:length(wl)
    x = sum(wl(1:w));
    plot([x x], [0 1000], 'color', 'w', 'linewidth', 3);
end

[W H] = nmfnnls(allcoh, 4);

subplot(3,4,9:12)
plot(H'); hold on
legend(lbl);
for w = 1:length(wl)
    x = sum(wl(1:w));
    plot([x x], [0 1000], 'color', 'k', 'linewidth', 3);
end
ylim([0 max(max(H))]);
xlim([1 size(H,2)]);

for p = 1:4
    subplot(3,4,4+p)
    imagesc(seeg_untril(W(:,p)));
    title(lbl{p});
    colorbar
    axis square
end

%% Identify optimal size of subgraphs into which to divide
%==========================================================================

[W H] = nnmf(dyncoh, 10);
subplot(3,1,1), imagesc(dyncoh);
subplot(3,1,2), imagesc(W*H);
subplot(3,1,3), imagesc(dyncoh - W*H);

textprogressbar('Simulating ');
count = 0;
for r = 1:10
for i = 1:50
    count = count + 1;
    textprogressbar(100*count/500);
    [W H]   = nnmf(dyncoh, i);
    err     = (dyncoh - W*H).^2;
    merr(r,i) = mean(mean(err));
end
end
textprogressbar(' Done');

figure
plot(mean(merr))
[val loc] = min(mean(merr));

%% 

opt = statset('maxiter', 100);
[W H] = nnmf(dyncoh, 6, 'algorithm', 'mult', 'replicates', 10, 'opt', opt);
opt = statset('maxiter', 1000, 'display', 'final');
[W H] = nnmf(dyncoh, 6, 'algorithm', 'als', 'w0', W, 'h0', H, 'opt', opt);
for i = 1:6
    wfac = max(W(:,i));
    W(:,i) = W(:,i) / wfac;
    H(i,:) = H(i,:) * wfac;
end


figure(1)
subplot(3,3,1:6)
    plot(H(1:3,:)');
    legend({'First', 'Second', 'Third'});

for p = 1:3
    subplot(3,3,6+p)
    imagesc(seeg_untril(W(:,p)));
    colorbar
    axis square
end

set(gcf, 'color', 'w')

figure(2)
subplot(3,3,1:6)
    plot(H(4:6,:)');
    legend({'First', 'Second', 'Third'});

for p = 1:3
    subplot(3,3,6+p)
    imagesc(seeg_untril(W(:,p+3)));
    colorbar
    axis square
end

set(gcf, 'color', 'w')


%% 
for i = 1:20
[W H it t r] = nmfnnls(dyncoh, i);
res(i)      = r;
end


