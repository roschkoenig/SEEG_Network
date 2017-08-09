% Housekeeping
%==========================================================================
clear all
D           = seeg_housekeeping;
Fanalysis   = D.Fanalysis;
chanlab     = D.chanlab;
fs          = filesep;
Fdcm        = [Fanalysis fs 'DCM'];
try mkdir(Fdcm); end
spm('defaults', 'EEG');
clear D


clear Ep
files = cellstr(spm_select('FPList', [Fdcm fs 'Individual Channels'], '.*.mat$'));
for f = 1:length(files)
    DCM = load(files{f});
    for i = 1:length(DCM.INV)
        INV = DCM.INV{i};
        Ep(f,i) = INV.Ep;
    end
end

% Time constants (one set)
%--------------------------------------------------------------------------
T   = mean(vertcat(Ep.T));

% Intrinsic connections (one per region)
%--------------------------------------------------------------------------
for f = 1:size(Ep,1)
    G(f,:) = mean(vertcat(Ep(f,:).G));
end

% Spatial model parameters (one set)
%--------------------------------------------------------------------------
L   = mean([Ep.L]);
J   = mean(vertcat(Ep.J));

% Region names
%--------------------------------------------------------------------------
for f = 1:length(files)
    seppos      = find(files{f} == fs);
    name{f}     = files{f}(seppos(end)+1:end-8);
end

P.T     = T;
P.G     = G;
P.L     = L;
P.J     = J;
P.name  = name;

save([Fdcm fs 'emp_priors'], 'P');