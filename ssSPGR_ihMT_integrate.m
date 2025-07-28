%%% [Mz_free] = ssSPGR_ihMT_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)
%
%   Steady-state ihMT SPGR sequence with eigenvector based time integration method
%
%   INPUTS:         
%           b1pulse    = RF pulse, Mx3 array (M=#timepoints, 3=frequency
%                        bands). Units are uT
%           dt         = dwell time, sec
%           Delta_Hz   = 1x3 vector of frequencies of each band. Typically 
%                        [-delta 0 delta]. Units Hz
%           TR         = repetition time, sec
%           tissuepars = structure containing all tissue parameters. See
%                        init_tissue()
%
%  OUTPUTS:
%           Mz_free         = Longitudinal magnetization of free water
%
% (c) Shaihan Malik 2019. King's College London

function [Mz_free] = ssSPGR_ihMT_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)

%%% Unpack tissue parameters
M0s = tissuepars.semi.M0;  %<--- TOTAL semisolid fraction
f   = tissuepars.semi.f;%<--- f = fraction of semisolid pool that is dipolar-coupled
M0f = 1-M0s;            %<--- free fraction is 1-M0s

R1  = [tissuepars.free.R1 tissuepars.semi.R1 tissuepars.semi.R1 tissuepars.semi.R1D];
R2f = tissuepars.free.R2;

% semisolid T2, used in lineshape calculation
T2s = tissuepars.semi.T2;

% Overall exchange rate for free and both semisolid pools
k = tissuepars.k;

%%% lineshape
switch tissuepars.lineshape
    case 'SL'
        [G,w_loc] = SuperLorentzian_lineshape(T2s,Delta_Hz,'interpzero');% seconds
    case 'Gaussian'
        [G,w_loc] = gauss_lineshape(T2s,Delta_Hz);% seconds
end

% gamma for RF calculation
gam = 267.5221; %< rad /s /uT



%% Lambda matrix and C are time invariant

% Free pool transverse components - no off resonance here
La = [-R2f 0;0 -R2f];

% the rest
Lb = [-k*M0s-R1(1) k*M0f k*M0f 0;k*(1-f)*M0s -k*M0f-R1(2) 0 0;...
    k*f*M0s 0 -k*M0f-R1(3) 0;0 0 0 -R1(4)];

Lambda = blkdiag(La,Lb);

C = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

% Augmented matrix form - this is a rate so applies whenever
Atilde = cat(1,[Lambda C],zeros(1,7));

% Work out pulse duration, and hence the duration of the evolution period.
% We can simply apply the constant rate for that whole period with no
% approximation
nt = size(b1pulse,1); %gives number of roles in b1pulse and shows how many discrete time-steps the pulse occupies
tau = dt*nt;


%% Now look at the RF matrices and integrate over pulse


% initialise this matrix
Xtilde_rf = eye(7);

for tt=1:nt
    
    % For simplicity assume that the pulses are all on x-axis 
    b1x = b1pulse(tt,2);% central band
    OmegaFree = gam*[0 0 0;0 0 b1x;0 -b1x 0];
    
    % now loop over the bands
    OmegaSemi = zeros(3,3);

    for ii=1:3
        w1 = gam*abs(b1pulse(tt,ii));
        w = pi*w1^2*G(ii);
        D = 2*pi*Delta_Hz(ii)/w_loc;
        OmegaSemi = OmegaSemi + [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]];
    end
    
    Omega = blkdiag(OmegaFree,OmegaSemi);
    
    % Now make overall evolution matrix
    Atilde = cat(1,[(Lambda+Omega) C],zeros(1,7));
    %% eigenvector decomposition or SVD, for eigenvalue = 0
    %% eigenvector is the steadystate
    
    % apply matrix exponential and multiply into existing matrix
    Xtilde_rf = expm(Atilde*dt)*Xtilde_rf;
    
end

%% steady-state with square pulse,get rid of the loop, do one run for one
%% small time increment
    

%% Now compile these and simulate longitudinal magnetization of free water

% initialize the M_0 vector with 7 elements
M0tilde = [ 0 ; 0 ; M0f ; (1-f)*M0s ; f*M0s ; 0 ; 1 ];
Mout = Xtilde_rf * M0tilde; % state right after the RF pulse
Mz_free = abs(Mout(3)); % the third element is the free-water longitudinal component 
        

end