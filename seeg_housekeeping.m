function D = seeg_housekeeping
% Housekeeping
%==========================================================================
fs          = filesep;

if strcmp(computer, 'PCWIN64') 
    Fbase = 'C:\Users\rrosch\Dropbox\Research\Friston Lab\1512 DCM collabs\sEEG DCM';
else
    Fbase = '/Users/roschkoenig/Dropbox/Research/Friston Lab/1512 DCM collabs/sEEG DCM'; 
end

Fscripts    = [Fbase fs 'Scripts'];
Fdata       = [Fbase fs 'Data'];
Fanalysis   = [Fbase fs 'Matlab Files'];
chanlab     = load([Fanalysis fs 'All_reduced_chanlist.mat']);
chanlab     = chanlab.reduced_channels;

spm('defaults', 'eeg');
addpath(genpath(Fscripts));

% Pack for exporting
%==========================================================================
D.Fbase     = Fbase;
D.Fscripts  = Fscripts;
D.Fdata     = Fdata;
D.Fanalysis = Fanalysis;
D.chanlab   = chanlab;
