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

%% Pre and post seizure PEB
%--------------------------------------------------------------------------
load([Fdcm fs 'DCM_Selection']);

si     = [3 4 5 6 7 8];
for s = 1:length(si)
    P{s} = INV{si(s)};
end
P = P';

Xnames  = {'Seizure', 'S1', 'S2', 'S3', 'GM'};
X       = [0 1 0 1 0 1;
           1 1 0 0 0 0; 
           0 0 1 1 0 0; 
           0 0 0 0 1 1;
           1 1 1 1 1 1]';

M.X         = X;
M.Xnames    = Xnames;

clear PEB
fields = {{'A'}, {'A', 'N'}, {'G'}, {'A','G'}};
for f = 1:length(fields)
    thisfield   = fields{f};
    PEB(f)      = spm_dcm_peb(P,M,thisfield);
end
[val loc]   = max([PEB.F]);

%% Run winning PEB and Bayesian Model Reduction
%==========================================================================
[PEB RCM]   = spm_dcm_peb(P,M,fields{loc});
BMA         = spm_dcm_peb_bmc(PEB);
save([Fdcm fs 'PMA_Seizures'], 'BMA');

%%
Ep      = BMA.Ep;
Cp      = diag(BMA.Cp);

ci      = spm_invNcdf(1 - 0.01);
c       = ci*sqrt(Cp);

Es      = abs(Ep) - c;          % High confidence parameters
Si      = find(Es > 0);

noP     = length(BMA.Pnames);
noC     = length(BMA.Xnames);

clear si ni
for c = 1:noC
    hi      = Si(Si >= (c-1)*noP);
    lo      = Si(Si <= c*noP);
    si{c}   = intersect(hi,lo);
    ni{c}   = PEB.Pnames(si{c} - (c-1)*noP);
end

%% Plot baseline network
%==========================================================================
A       = zeros(7);
ci      = spm_invNcdf(1 - 0.01);
c       = ci*sqrt(Cp);
node    = zeros(7,1);

% Collate significant edges
%--------------------------------------------------------------------------
for p = 1:length(BMA.Pnames)
    N = BMA.Pnames{p};

    TEp = BMA.Ep(p + noP*4);
    Tci = c(p+noP*4);

    if abs(TEp) - Tci > 0    

        % Identify significant extrinsic connection changes
        %------------------------------------------------------------------
        if N(1) == 'A'
            ri = str2double(N(6));
            ci = str2double(N(8));

            A(ri,ci)    = A(ri,ci) + TEp;

        % Identify significant intrinsic connection changes
        %------------------------------------------------------------------
        elseif N(1) == 'G'
            nodeid          = str2double(N(3));
            node(nodeid)    = TEp;
        end
    end
end

cols        = flip(cbrewer('div', 'Spectral', 100));
colormap(cols)

% Make directed graph object
%--------------------------------------------------------------------------
G           = digraph(A, RCM{1}.Sname);

% Plot directed graph
%--------------------------------------------------------------------------
LWidths     = abs(5*G.Edges.Weight/max(G.Edges.Weight));
ECdata      = 10*G.Edges.Weight ./ abs(G.Edges.Weight);
NCdata      = node;
plot(G, 'LineWidth', LWidths, 'Layout', 'circle', 'EdgeCData', ECdata, 'NodeCData', NCdata)
caxis([-1 1]);
axis square


%% Plot network connectivity changes for seizure conditions
%==========================================================================
for cond = 1:3
    condid = cond + 1;
    
    A       = zeros(7);
    ci      = spm_invNcdf(1 - 0.01);
    c       = ci*sqrt(Cp);

    % Collate significant edges
    %----------------------------------------------------------------------
    node = zeros(7,1);
    
    for p = 1:length(BMA.Pnames)
        N = BMA.Pnames{p};

        TEp = BMA.Ep(p + noP*(condid-1));
        Tci = c(p+noP*(condid-1));

        if abs(TEp) - Tci > 0    

            % Identify significant extrinsic connection changes
            %------------------------------------------------------------------
            if N(1) == 'A'
                ri = str2double(N(6));
                ci = str2double(N(8));

                A(ri,ci)    = A(ri,ci) + TEp;

            % Identify significant intrinsic connection changes
            %------------------------------------------------------------------
            elseif N(1) == 'G'
                nodeid          = str2double(N(3));
                node(nodeid)    = TEp;
            end
        end
    end

    cols    = flip(cbrewer('div', 'Spectral', 100));
    colormap(cols)
    
    % Make directed graph obeject
    %----------------------------------------------------------------------
    G           = digraph(A, RCM{1}.Sname);

    % Plot directed graph
    %----------------------------------------------------------------------
    LWidths     = abs(5*G.Edges.Weight/max(G.Edges.Weight));
    ECdata      = G.Edges.Weight ./ abs(G.Edges.Weight);
    NCdata      = node;
    subplot(1,3,cond)
    plot(G, 'LineWidth', LWidths, 'Layout', 'circle', 'EdgeCData', ECdata, 'NodeCData', NCdata)
    caxis([-1 1]);
    axis square
    
end

%% Plot intrinsic connectivity changes leading to seizure
%==========================================================================
condid = 1;

A       = zeros(7);
ci      = spm_invNcdf(1 - 0.01);
c       = ci*sqrt(Cp);

% Collate significant edges
%----------------------------------------------------------------------
node = zeros(7,1);

for p = 1:length(BMA.Pnames)
    N = BMA.Pnames{p};

    TEp = BMA.Ep(p + noP*(condid-1));
    Tci = c(p+noP*(condid-1));

    if abs(TEp) - Tci > 0    

        % Identify significant extrinsic connection changes
        %------------------------------------------------------------------
        if N(1) == 'A'
            ri = str2double(N(6));
            ci = str2double(N(8));

            A(ri,ci)    = A(ri,ci) + TEp;
        
        % Identify significant intrinsic connection changes
        %------------------------------------------------------------------
        elseif N(1) == 'G'
            nodeid          = str2double(N(3));
            node(nodeid)    = TEp;
        end
    end
end

cols    = flip(cbrewer('div', 'Spectral', 100));
colormap(cols)

% Make directed graph obeject
%----------------------------------------------------------------------
G           = digraph(A, RCM{1}.Sname);

% Plot directed graph
%----------------------------------------------------------------------
LWidths     = abs(5*G.Edges.Weight/max(G.Edges.Weight));
ECdata      = G.Edges.Weight ./ abs(G.Edges.Weight);
NCdata      = node;

plot(G, 'LineWidth', LWidths, 'Layout', 'circle', 'EdgeCData', ECdata, 'NodeCData', NCdata);
caxis([-1 1]);
axis square






