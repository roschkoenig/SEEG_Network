D         = seeg_housekeeping;
Fbase     = D.Fbase;
Fscripts  = D.Fscripts;
Fdata     = D.Fdata;
Fanalysis = D.Fanalysis;
fs = filesep;

for s = 1:length(S)
    subplot(5,1,s)
    plot(S(s).dat(10,:));
    title(S(s).name);
    set(gcf, 'Color', 'w');
end
%%
C = S(2);
win     = 30 * C.hdr.Fs;
step    = fix(0.1 * win);

%% Calculate sliding window correlation and bandpower distribution
%--------------------------------------------------------------------------
clear count thiswin wdat wcor wpow wcoh
count   = 0;
textprogressbar('Sliding Window: ');
for w = 1 : step : C.hdr.nSamples - win
    count    = count + 1;
    thiswin  = [0 : win - 1] + w;
    wdat     = C.dat(:, thiswin);
    wcor(count,:,:)  = corr(wdat');
    wcoh(count,:,:)  = seeg_coh(wdat, 10);
    for c = 1:size(C.dat,1)
       wpow(count,c) = bandpower(wdat(c,:));
    end
	textprogressbar(count * 100 / length(1 : step : C.hdr.nSamples - win))
end
textprogressbar(' Done');

% Save estimates in a structure
%==========================================================================
B.wpow = wpow;
B.wcor = wcor;
B.wcoh = wcoh;

save([Fanalysis fs C.name '_win'], 'B');

%% Plot bandpower distribution over time, and network correlation 
%==========================================================================
figure(1)
subplot(5,1,[1:4]);
    imagesc(wpow'), colorbar; colormap parula;
    set(gca, 'YTick', 1:size(C.dat,1));
    set(gca, 'YTickLabel', C.hdr.label(2:size(C.dat,1)+1));
subplot(5,1,5)
    mfg08 = find(strcmp('MFG08', C.hdr.label)) - 1;
    plot(C.dat(mfg08,:));
set(gcf, 'Color', 'w');

figure(2)
subplot(2,2,1)
    imagesc(squeeze(wcor(30,:,:))); axis square
    title('Correlation baseline');
	set(gca, 'YTick', 1:size(C.dat,1));
    set(gca, 'YTickLabel', C.hdr.label(2:size(C.dat,1)+1));
subplot(2,2,2)
    imagesc(squeeze(wcor(150,:,:))); axis square
    title('Correlation seizure');
	set(gca, 'YTick', 1:size(C.dat,1));
    set(gca, 'YTickLabel', C.hdr.label(2:size(C.dat,1)+1));
    
subplot(2,2,3)
    imagesc(squeeze(wcoh(30,:,:))); axis square
    title('Coherence baseline');
	set(gca, 'YTick', 1:size(C.dat,1));
    set(gca, 'YTickLabel', C.hdr.label(2:size(C.dat,1)+1));
subplot(2,2,4)
    imagesc(squeeze(wcoh(150,:,:))); axis square
    title('Coherence seizure');
	set(gca, 'YTick', 1:size(C.dat,1));
    set(gca, 'YTickLabel', C.hdr.label(2:size(C.dat,1)+1));
set(gcf, 'Color', 'w');

%% Plot corrrelation and coherence changes over time
%==========================================================================
for c = 1:size(wcor,1)
    thiscor      = squeeze(wcor(c,:,:));
    triid        = find(tril(ones(size(wcor,2)),-1));     % lower triangle indices
    dyncor(:,c)  = thiscor(triid);
end
mcor    = mean(dyncor,2);
[sd si] = sort(mcor);
subplot(2,1,1)
imagesc(dyncor(si,:));
set(gcf, 'color', 'w');

for c = 1:size(wcoh,1)
    thiscoh      = squeeze(wcoh(c,:,:));
    triid        = find(tril(ones(size(wcoh,2)),-1));     % lower triangle indices
    dyncoh(:,c)  = thiscoh(triid);
end
mcoh    = mean(dyncoh,2);
[sd si] = sort(mcoh);
subplot(2,1,2)
imagesc(log(dyncoh(si,:)));
set(gcf, 'color', 'w');

% Plot dynamics matrices
%==========================================================================
cpow = corr(wpow'); subplot(1,3,1), imagesc(cpow); axis square; title('Power')
ccor = corr(dyncor); subplot(1,3,2), imagesc(ccor); axis square; title('Correlation');
ccoh = corr(dyncoh); subplot(1,3,3), imagesc(ccoh); axis square; title('Coherence');
set(gcf, 'color', 'w');

