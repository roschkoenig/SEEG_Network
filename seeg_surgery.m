% SEEG Surgical Simulation
%==========================================================================
% This routine uses the inverted DCMs of the peri-seizure windows from the
% previous analysis, and simulates surgical interention based on those to
% predict the impact of seizures on both coherence and spectral densities 

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

% Pre and post seizure PEB
%==========================================================================
load([Fdcm fs 'DCM_Selection']);

% Collate DCMs into single matrix
%--------------------------------------------------------------------------
si     = [3 4 5 6 7 8];
for s = 1:length(si)
    P{s} = INV{si(s)};
end
P = P';

% Set up PEB regressors
%--------------------------------------------------------------------------
Xnames  = {'Seizure', 'S1', 'S2', 'S3', 'GM'};
X       = [0 1 0 1 0 1;
           1 1 0 0 0 0; 
           0 0 1 1 0 0; 
           0 0 0 0 1 1;
           1 1 1 1 1 1]';

M.X         = X;
M.Xnames    = Xnames;

% Run PEB and derive reduced models
%--------------------------------------------------------------------------
clear PEB
[PEB RCM]   = spm_dcm_peb(P,M,{'A', 'G'});

% Loop through removal of different nodes
%==========================================================================
for s = 1:7
surgnode = s;

for ri = [2 4 6]
    int_Ep      = RCM{ri}.Ep;
    surg_Ep     = int_Ep;
    A           = int_Ep.A;
    
    % Disconnect the surgical node
    %----------------------------------------------------------------------
    Ns  = length(A{1});
    for a = 1:length(A)
        A{a}(surgnode,:)    = ones(length(surgnode),Ns)*(-32);
        A{a}(:,surgnode)    = ones(Ns,length(surgnode))*(-32);
    end
    surg_Ep.A   = A;
    
    % Estimate the resultant cross-spectral densities
    %----------------------------------------------------------------------
    intact_csd  = spm_csd_mtf(int_Ep, RCM{ri}.M, RCM{ri}.xU);
    surg_csd    = spm_csd_mtf(surg_Ep, RCM{ri}.M, RCM{ri}.xU); 
    
    % Calculate difference in CSD and Coherence
    %----------------------------------------------------------------------
    diff_csd    = abs(intact_csd{1}) - abs(surg_csd{1});
    diff_coh    = seeg_csd2coh(intact_csd{1}) - seeg_csd2coh(surg_csd{1});
    mdiffcsd(s)    = mean(real(diff_csd(:)));
    mdiffcoh(s)    = mean(diff_coh(:));
end
end

% Plot results
%--------------------------------------------------------------------------
cols = [227,6,19; 249,175,21; 48,146,54; 0,144,199; 56,170,57; 83,55,141; 120,68,148] / 255;
mdiffcsd = abs(mdiffcsd) / max(abs(mdiffcsd));
mdiffcoh = abs(mdiffcoh) / max(abs(mdiffcoh));


Ns      = length(mdiffcsd);
xpos    = ones(1,Ns) + (rand(1,Ns)-.5)/3 ;
scatter(xpos, mdiffcoh, 150, cols, 'filled'); hold on

xpos    = ones(1,Ns) + (rand(1,Ns)-.5)/3 + 1;
scatter(xpos, mdiffcsd, 150, cols, 'filled');

xlim([0 3]);
plot([1 1], [0 1], 'color', 'k');
plot([2 2], [0 1], 'color', 'k');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Spectral difference', 'Coherence Difference'});
ylabel('Normalised difference induced by surgery')
set(gcf, 'color', 'w');
