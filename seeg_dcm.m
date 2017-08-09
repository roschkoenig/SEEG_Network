% function [INV] = seeg_dcm(sub, Fanalysis, invert) 

% Housekeeping
%==========================================================================
clear DCM
D           = seeg_housekeeping;
Fanalysis   = D.Fanalysis;
chanlab     = D.chanlab;
fs          = filesep;
Fdcm        = [Fanalysis fs 'DCM'];
try mkdir(Fdcm); end

spm('defaults', 'EEG');

% Specify analysis segment
%==========================================================================
win_no = [150 1105 323 343 710 730 865 885];

% Setup DCM Structure
%==========================================================================
DCM = [];

% Fix directory of canonical forward matrix
%--------------------------------------------------------------------------
DCM.xY.Dfile      	= [Fanalysis fs 'MEEG' fs 'SEEG.mat'];

% Load MEEG object and extract sampling rate and info
%--------------------------------------------------------------------------
SEEG                = spm_eeg_load(DCM.xY.Dfile);
Fs                  = fsample(SEEG);
smpls               = size(SEEG,2);
timax               = linspace(0, smpls/Fs, smpls);
conds               = condlist(SEEG);

for w = 1:length(win_no)
    thiscond            = conds{win_no(w)};
    display(thiscond)

    % Set up DCM details
    %--------------------------------------------------------------------------
    DCM.options.analysis = 'CSD';       % cross-spectral density 
    DCM.options.model    = 'CMC';      	% structure canonical microcircuit (for now)
    DCM.options.spatial  = 'LFP';           
    DCM.options.Tdcm     = [timax(1) timax(end)] * 1000;     

    DCM.options.Fdcm    = [1 60];     	% frequency range  
    DCM.options.D       = 1;         	% frequency bin, 1 = no downsampling
    DCM.options.Nmodes  = 8;          	% cosine reduction components used 
    DCM.options.han     = 0;         	% no hanning 
    DCM.options.trials  = w;            % index of ERPs within file

    DCM.Sname           = chanlabels(SEEG);
    DCM.M.Hz            = DCM.options.Fdcm(1):DCM.options.D:DCM.options.Fdcm(2);
    DCM.xY.Hz           = DCM.M.Hz;

    % Define different connection types
    %==========================================================================
    F      = tril(ones(size(SEEG,1)),-1);    % Forward triangle (lower)
    B      = triu(ones(size(SEEG,1)),1);     % Backward triangle (upper)
    S      = zeros(size(SEEG,1));   for s = 1:length(S); S(s,s) = 1; end

    % Define model arcitecture (A), conditional effects (B) and input (C) 
    %==========================================================================
    DCM.A{1}    =   F + B;
    DCM.A{2}    =   F + B;
    DCM.A{3}    =   S;

    DCM.B           = {};
    DCM.C           = sparse(length(DCM.A{1}),0); 

    % Reorganise model parameters in specific structure
    %==========================================================================
    DCM.M.dipfit.Nm    = DCM.options.Nmodes;
    DCM.M.dipfit.model = DCM.options.model;
    DCM.M.dipfit.type  = DCM.options.spatial;

    DCM.M.dipfit.Nc    = size(SEEG,1);
    DCM.M.dipfit.Ns    = length(DCM.A{1});

    % Define priors
    %==========================================================================
    % Load standard neural priors
    %--------------------------------------------------------------------------
    [pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);  
    [pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
    [pE,pC]  = spm_ssr_priors(pE,pC);
    
    load([Fdcm fs 'emp_priors.mat']);
    pE.T        = P.T;
    pE.L        = ones(length(pE.L),1) .* P.L;
    pE.J        = P.J;

    for n = 1:length(P.name)
        thisname        = P.name{n};
        nameid          = find(strcmp(chanlab, thisname));
        pE.G(nameid,:)  = P.G(n,:);
    end

    % Switch off variations in spatial parameters
    %--------------------------------------------------------------------------
    pC.L        = pC.L * 0;
    pC.J        = pC.J * 0;

    DCM.M.pE   = pE;
    DCM.M.pC   = pC;


    % Invert DCM and log inversion output
    %==========================================================================
    diary([Fdcm fs 'log']);
    DCM.name    = [Fdcm fs 'DCM_' conds{win_no(w)}];
    TMP         = seeg_spm_dcm_csd(DCM);
    TMP.xY.R    = diag(TMP.xY.R);
    INV{w}      = TMP;
    
	delete(DCM.name);
    save([Fdcm fs 'DCM_Selection'], 'INV');
    diary off
end
