clear all
D         = seeg_housekeeping;
Fbase     = D.Fbase;
Fscripts  = D.Fscripts;
Fdata     = D.Fdata;
Fanalysis = D.Fanalysis;
fs        = filesep;

%% Load data into single structure
%==========================================================================
files = cellstr(spm_select('FPList', Fdata, '^*.edf'));
for f = 1:length(files)
    file = files{f};
    seps        = find(file == fs);
    S(f).name   = file(seps(end)+1:end-4);
    S(f).hdr    = ft_read_header(file);
    S(f).dat    = ft_read_data(file);
    chid        = find(strcmp('PAR08', S(f).hdr.label));
    S(f).dat    = S(f).dat(2:chid,:);
end

%% Save data as MEEG object
%==========================================================================
count   = 0;
for s = 1:length(S)
C = S(s);

load([Fanalysis fs 'All_reduced_chanlist']);

for c = 1:length(reduced_channels)
    ch      = reduced_channels{c};
    chi(c)  = find(strcmp(ch, C.hdr.label));
end

nS      = C.hdr.nSamples;
Fs      = C.hdr.Fs;
win     = 30 * Fs;
step    = fix(0.1 * win);

for w = 1:step:(nS-win)
    count = count + 1;
    ftdata.trial{count} = C.dat(chi, w:w+win);
    ftdata.time{count}  = [1:1+win] / Fs;
    conds{count}        = [C.name '_' num2str(fix(w / Fs))];
end
end
ftdata.label = C.hdr.label(chi);

D = spm_eeg_ft2spm(ftdata, [Fanalysis fs 'MEEG' fs 'SEEG']);

for d = 1:size(D,3)
    D = conditions(D, d, conds{d});
end
save(D);

%% MEEG with single channel conditions
%==========================================================================
count   = 0;
data    = [];
ftdata  = [];

% Loop through individual files
%--------------------------------------------------------------------------
for s = 1:length(S)
C = S(s);

load([Fanalysis fs 'All_reduced_chanlist']);

for c = 1:length(reduced_channels)
    ch      = reduced_channels{c};
    chi(c)  = find(strcmp(ch, C.hdr.label));
end

nS      = C.hdr.nSamples;
Fs      = C.hdr.Fs;

data    = [data, C.dat(chi,:)];

end

% Split into individual channels
%--------------------------------------------------------------------------
for d = 1:size(data,1)
    ftdata.trial{d} = data(d,:);
    ftdata.time{d}  = [1:size(data,2)] / Fs;
end
ftdata.label    = {'Channel'};
D               = spm_eeg_ft2spm(ftdata, [Fanalysis fs 'MEEG' fs 'SEEG_cont_by_chan']);

for d = 1:size(D,3)
    D = conditions(D, d, C.hdr.label{chi(d)});
end
save(D);
