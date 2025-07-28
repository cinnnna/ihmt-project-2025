function [Mz_free] = ssSPGR_ihMT_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)
%%% [Mz_free] = ssSPGR_ihMT_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)
%
% Steady-state ihMT SPGR sequence with eigenvector method for square pulses
%
% INPUTS:
% b1pulse = RF pulse, Mx3 array (M=#timepoints, 3=frequency
% bands). Units are uT. For square pulse, all rows are identical.
% dt = dwell time, sec
% Delta_Hz = 1x3 vector of frequencies of each band. Typically
% [-delta 0 delta]. Units Hz
% tissuepars = structure containing all tissue parameters. See
% init_tissue()
%
% OUTPUTS:
% Mz_free = Longitudinal magnetization of free water
%
% Modified for square pulse using eigenvector method

%%% Unpack tissue parameters
M0s = tissuepars.semi.M0; %<--- TOTAL semisolid fraction
f = tissuepars.semi.f;%<--- f = fraction of semisolid pool that is dipolar-coupled
M0f = 1-M0s; %<--- free fraction is 1-M0s
R1 = [tissuepars.free.R1 tissuepars.semi.R1 tissuepars.semi.R1 tissuepars.semi.R1D];
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

%% For square pulse, RF is constant - calculate Omega once
% Take first row of b1pulse (all rows should be identical for square pulse)
b1x = b1pulse(1,2);% central band
OmegaFree = gam*[0 0 0;0 0 b1x;0 -b1x 0];

% now loop over the bands for semisolid pool
OmegaSemi = zeros(3,3);
for ii=1:3
    w1 = gam*abs(b1pulse(1,ii));
    w = pi*w1^2*G(ii);
    D = 2*pi*Delta_Hz(ii)/w_loc;
    OmegaSemi = OmegaSemi + [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]];
end

Omega = blkdiag(OmegaFree,OmegaSemi);

%% Form the augmented matrix
% Augmented matrix form for homogeneous equation
Atilde = cat(1,[(Lambda+Omega) C],zeros(1,7));

%% Calculate steady-state using eigenvector method
% Find eigenvector corresponding to eigenvalue 0
[V,D] = eig(Atilde);

% Find the eigenvalue with minimum absolute value
eigenvals = diag(D);
[min_eval,ind] = min(abs(eigenvals));

% Verify it's actually zero (within machine precision)
if min_eval > eps*100
    error('No zero eigenvalue found. Minimum eigenvalue is %e. Check system setup.', min_eval);
end

% Extract the corresponding eigenvector
Mss = V(:,ind);

% Normalize by the augmented element (last element should be 1)
if abs(Mss(end)) > eps
    Mss = Mss/Mss(end);
else
    warning('Augmented element is too small, using alternative normalization');
    % Alternative: normalize by total magnetization
    Mss = Mss/sum(abs(Mss(3:5)));
end

% Extract free water longitudinal magnetization
Mz_free = real(Mss(3));

% Ensure result is physical (between 0 and M0f)
Mz_free_raw = max(0, min(M0f, Mz_free));
Mz_free = Mz_free_raw / M0f; %Mz normalized to M0

end