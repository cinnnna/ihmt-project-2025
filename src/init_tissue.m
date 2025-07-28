%%% pars = init_tissue(name)
%
% Function to generate structure containing tissue properties. Possible
% name inputs are:
%
% - 'ic' = 'Internal Capsule'
% - 'simple' = simple rounded parameter values
% - 'pl161' = prolipid 161 as measured experimentally
%
% default = 'ic'
%
% The function initialises the following parameters:
%%% Feb 2019 Initialise tissue parameters
function pars = init_tissue(name)
pars = struct;
if nargin==0
    % set to default
    name = 'ic';
end
% Replace all letters with lower case
name = lower(name);
% initialise the empty structure
% Rates are all s^-1, times on seconds, magnetizations referenced to total
% M0=1
%%% Free water pool
pars.free.R1 = [];
pars.free.R2 = [];
%%% semisolid pool
pars.semi.M0 = []; % assume M0f+M0s = 1
pars.semi.R1 = []; %<--- assume same for both s1 and s2
pars.semi.R1D = [];
pars.semi.f = []; % fraction of M0s that is in s1
pars.semi.T2 = []; %<--- assume same for both s1 and s2
%%% overall exchange rate constant
pars.k = [];

switch name
    case 'random'
        rng(42); % or use rng('shuffle') for different results each time
        % Define reasonable ranges for each parameter based on literature
        ranges = struct();
        ranges.free_R1 = [0.1, 2.0];      % 1/T1 for free water (s^-1)
        ranges.free_R2 = [1, 50];         % 1/T2 for free water (s^-1) 
        ranges.semi_M0 = [0.01, 0.3];     % Semisolid fraction
        ranges.semi_R1 = [1, 20];         % Semisolid R1 (s^-1)
        ranges.semi_R1D = [10, 100];      % Dipolar R1 (s^-1)
        ranges.semi_T2 = [5e-6, 50e-6];   % Semisolid T2 (s)
        ranges.k = [1, 20];             % Exchange rate (s^-1)

        % Generate random initial values
        pars.free.R1 = ranges.free_R1(1) + diff(ranges.free_R1) * rand();
        pars.free.R2 = ranges.free_R2(1) + diff(ranges.free_R2) * rand();
        pars.semi.M0 = ranges.semi_M0(1) + diff(ranges.semi_M0) * rand();
        pars.semi.R1 = ranges.semi_R1(1) + diff(ranges.semi_R1) * rand();
        pars.semi.R1D = ranges.semi_R1D(1) + diff(ranges.semi_R1D) * rand();
        pars.semi.f = 0; % Keep this fixed %1 hc 0 eggwhite
        pars.semi.T2 = ranges.semi_T2(1) + diff(ranges.semi_T2) * rand();
        pars.k = ranges.k(1) + diff(ranges.k) * rand();
        pars.lineshape = 'Gaussian';

    case 'bsa'
        pars.free.R1 = 0.5635;
        pars.free.R2 = 18.812;
        pars.semi.M0 = 0.0774;
        pars.semi.R1 = 7.7545;
        pars.semi.R1D = 47.8069;
        pars.semi.f = 0;
        pars.semi.T2 = 15.20e-6;
        %%% overall exchange constant
        pars.k = 70.1;
        pars.lineshape = 'Gaussian';
        
    case 'hc'
        pars.free.R1 = 0.5955; % s⁻¹ (T1f ≈ 1.68 s)  %random value
        pars.free.R2 = 5.9999; % s⁻¹ (T2f ≈ 167 ms)
        % ---- Semisolid (Zeeman + dipolar) pools -----------------------------
        pars.semi.M0 = 0.0675; % total semisolid fraction (1–M0f)
        pars.semi.R1 = 2.9780; % s⁻¹ (Zeeman pool)
        pars.semi.R1D = 42.8042; % s⁻¹ (dipolar-order pool)
        pars.semi.f = 0.7139; % fraction that is dipolar-coupled
        pars.semi.T2 = 25.85e-6; % s (25.9 µs)
        pars.k = 67.8720; % exchange rate
        % ---- Lineshape ------------------------------------------------------
        pars.lineshape = 'SL'; % Super-Lorentzian
        
    case 'simple'
        pars.free.R1 = 1;
        pars.free.R2 = 10; % 100ms T2
        pars.semi.M0 = 0.2;
        pars.semi.R1 = 1;
        pars.semi.R1D = 1/10e-3;
        pars.semi.f = 0.65;
        pars.semi.T2 = 12e-6;
        pars.k = 50;
        pars.lineshape = 'SL'; %<-- super Lorentzian
        
    case 'ic'
        %%% Parameters from Mchinda 2017 -- Default parameters
        pars.free.R1 = 1/650e-3;
        pars.free.R2 = 1/80e-3;
        %pars.semi.M0 = 0.17; %<--- assume that total M0 is 1
        pars.semi.M0 = 0.147; %<--- M0f=1 for Varma work. Here M0s = M0s_Varma / (1+M0s_Varma)
        pars.semi.R1 = 1; %<--- assume same for both s1 and s2
        pars.semi.R1D = 1/6.5e-3; %<--- dipolar relaxation rate
        pars.semi.f = 0.65; %<--- fraction that has long T1D
        pars.semi.T2 = 12e-6; %<--- assume same for both s1 and s2
        %%% overall exchange constant
        pars.k = 65;
        pars.lineshape = 'SL';
        
    case 'pl161' %<-- fitted data from 2019 paper
        pars.free.R1 = 1/1764e-3;
        pars.free.R2 = 1/73e-3;
        pars.semi.M0 = 0.151;
        pars.semi.R1 = 1/226e-3; %<--- assume same for both s1 and s2
        pars.semi.R1D = 1/20.7e-3; %<--- dipolar relaxation rate
        pars.semi.f = 1; %<--- fraction that has long T1D
        pars.semi.T2 = 17.9e-6; %<--- assume same for both s1 and s2
        %%% overall exchange constant
        pars.k = 46.7;
        pars.lineshape = 'Gaussian';
        
    otherwise
        error('Unknown tissue type: %s', name);
end

end