function [Mz_free] = Dualcase_ssSPGR_ihMT_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)
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
% Take the on-resonance component for free pool
% For free pool, we need the pulse amplitude at 0 Hz (on-resonance)
zero_idx = find(abs(Delta_Hz) < eps); % Find which band is at 0 Hz
if ~isempty(zero_idx)
    b1x = b1pulse(1, zero_idx(1)); % Use pulse from the 0 Hz band
else
    b1x = 0; % No on-resonance component
end

OmegaFree = gam*[0 0 0;0 0 b1x;0 -b1x 0];

% now loop over the bands for semisolid pool - 2-BAND SYSTEM
OmegaSemi = zeros(3,3);

% Always process 2 bands: [-delta, 0], [-delta, +delta], or [0, +delta]
for ii = 1:2
    w1 = gam*abs(b1pulse(1,ii));
    w = pi*w1^2*G(ii);
    
    % Special case: if Delta_Hz(ii) = 0, then D = 0 (w*D = 0)
    if abs(Delta_Hz(ii)) < eps
        D = 0; % w*D = 0 when delta = 0
    else
        D = 2*pi*Delta_Hz(ii)/w_loc;
    end
    
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