% Identifying subset of channels from whole montage
%==========================================================================
% this routine uses nonnegative matrix decomposition to identify a subset
% of channels correspond most to highly coherent subsections of the montage

% Housekeeping
%==========================================================================
clear all
D         = seeg_housekeeping;
Fbase     = D.Fbase;
Fscripts  = D.Fscripts;
Fdata     = D.Fdata;
Fanalysis = D.Fanalysis;
fs = filesep;

% Load all available data into single structure
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

% Perform nonnegative matrix decomposition
%==========================================================================
[w h] = nnmf(allcoh, 5);

%% Plot decomposition results
%--------------------------------------------------------------------------
subplot(2,1,1); plot(w)
subplot(2,1,2); plot(h')
set(gca, 'xtick', 1:38);
set(gca, 'xticklabel', chanlab);
figure(1)
subplot(2,2,4), 
    imagesc(allcoh), axis square;
    title('Original Matrix')
subplot(2,2,3),
    imagesc(w);
    ylabel('channels'); xlabel('components');
	title('W Components');
subplot(2,2,2),
    imagesc(h);
    title('H Components');
    xlabel('channels'); ylabel('components');
    set(gcf, 'color', 'w');
    xlabel('channels'); ylabel('components');

% Identify channels containing most representative of local coherence
%--------------------------------------------------------------------------
mw = max(w');
mh = max(h);
figure(2)
subplot(2,1,1), 
    plot(w); hold on, 
    plot(mw, 'r');
    title('Expression of W Components');
    legend({'Component 1', 'Component 2', 'Component 3', 'Component 4', 'Component 5', 'Envelope'});
    ylabel('Expression');
    xlabel('Channel');
    set(gcf, 'color', 'w');
    
subplot(2,1,2), 
    findpeaks(mw);
    title('Peaks of Expression');

    [vw lw] = findpeaks(mw);
    [vh lh] = findpeaks(mh);
    
    set(gca, 'XTick', lw);
    set(gca, 'XTickLabel', chanlab(lw)); 

% Save identified channel labels
%==========================================================================
channels = unique([lw lh]);
reduced_channels = chanlab(channels);
save([Fanalysis fs 'All_reduced_chanlist'], 'reduced_channels');


