function [Mz_free] = new_Dualcase_ssSPGR_ihMT_integrate(b1_band_value,dt,delta,tissuepars,varargin)
%%% Modified version that normalizes like experimental data
%%% Divides by M0 calculated with B1_max = 0 (no RF pulse)

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
gam = 267.5221; %< rad /s /uT

%% Lambda matrix and C are time invariant
% Free pool transverse components - no off resonance here
La = [-R2f 0;0 -R2f];
% the rest
Lb = [-k*M0s-R1(1) k*M0f k*M0f 0;k*(1-f)*M0s -k*M0f-R1(2) 0 0;...
    k*f*M0s 0 -k*M0f-R1(3) 0;0 0 0 -R1(4)];
Lambda = blkdiag(La,Lb);
C = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

%% Calculate Omega matrices
b1x = real(b1_band_value); % Real part (x-component) of complex pulse
b1y = imag(b1_band_value); % Imaginary part (y-component) of complex pulse

% fprintf('b1x = %.6f, b1y = %.6f\n', b1x, b1y);
% fprintf('Pulse shape (first value): %.3f + %.3fi\n', real(b1_band(1)), imag(b1_band(1)));
OmegaFree = gam*[0 0 -b1y;0 0 b1x;b1y -b1x 0];

% Semisolid pool
    %%% lineshape
    switch tissuepars.lineshape
        case 'SL'
            [G,w_loc] = SuperLorentzian_lineshape(T2s,delta,'interpzero');% seconds
        case 'Gaussian'
            [G,w_loc] = gauss_lineshape(T2s,delta);% seconds
    end

% Process single complex band
w1 = gam*abs(b1_band_value);  % RF amplitude from complex pulse magnitude
w = pi*w1^2*G;          % Use first element of lineshape (assuming single frequency)


% Calculate frequency offset effect
if abs(delta) < eps
    D = 0; % w*D = 0 when on-resonance
else
    D = 2*pi*delta/w_loc;
end


OmegaSemi = [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]];

Omega = blkdiag(OmegaFree,OmegaSemi);

%% Form the augmented matrix
% Augmented matrix form for homogeneous equation
Atilde = cat(1,[(Lambda + Omega) C],zeros(1,7));

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


% Normalize the engenvector so that the last component is 1
%Mss = Mss / Mss(end);

% Extract free water longitudinal magnetization (RAW value, not normalized to M0f)
Mz_free = real(Mss(3))/M0f;

end