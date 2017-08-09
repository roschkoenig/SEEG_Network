clear all
D         = seeg_housekeeping;
Fbase     = D.Fbase;
Fscripts  = D.Fscripts;
Fdata     = D.Fdata;
Fanalysis = D.Fanalysis;
fs = filesep;

%% Load data into single structure
%==========================================================================
files = cellstr(spm_select('FPList', Fdata, '^*.edf'));
for f = 1:length(files)
    file = files{f};
    S(f).name   = file(77:end-4);
    S(f).hdr    = ft_read_header(file);
    S(f).dat    = ft_read_data(file);
    chid        = find(strcmp('PAR08', S(f).hdr.label));
    S(f).dat    = S(f).dat(2:chid,:);
end

chanlab = S(1).hdr.label(2:39);
alldat  = [];

for s = 1:length(S)
    alldat = [alldat, S(s).dat];
end
allcoh = seeg_coh(alldat,20);
%%

[w h] = nnmf(allcoh, 4);
subplot(2,1,1); plot(w)
subplot(2,1,2); plot(h')
set(gca, 'xtick', 1:38);
set(gca, 'xticklabel', chanlab);

subplot(2,2,4), imagesc(allcoh), axis square;
subplot(2,2,3), imagesc(w);
subplot(2,2,2), imagesc(h);

mw = max(w');
mh = max(h);
subplot(2,1,1), plot(w); hold on, plot(mw, 'r');
subplot(2,1,2), findpeaks(mw);
[vw lw] = findpeaks(mw);
[vh lh] = findpeaks(mh);

channels = unique([lw lh]);
reduced_channels = chanlab(channels);
save([Fanalysis fs 'All_reduced_chanlist'], 'reduced_channels');


