function [Mz_free] = new_Dualcase_ssSPGR_ihMT_integrate(b1_band,dt,delta,tissuepars,nband)
%%% Modified version using b1_band and delta
%%% Integrates over all entries in b1_band vector

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

%% FIRST: Calculate M0 with no RF pulse (B1 = 0)
% Create zero pulse matrix
b1pulse_zero = zeros(size(b1_band));
Delta_Hz_zero = 0; % No frequency offset matters when B1=0

% Calculate reference state (no RF)
Mz_0 = calculate_steady_state(b1pulse_zero, dt, Delta_Hz_zero, tissuepars,nband, M0s, f, M0f, R1, R2f, T2s, k);

%% SECOND: Calculate Mz with actual RF pulse
Mz_with_RF = calculate_steady_state(b1_band, dt, delta, tissuepars,nband, M0s, f, M0f, R1, R2f, T2s, k);

%% THIRD: Normalize like experimental data (Mz/M0)
Mz_free = Mz_with_RF / Mz_0;

end

function Mz_out = calculate_steady_state(b1_band, dt, delta, tissuepars,nband, M0s, f, M0f, R1, R2f, T2s, k)
% gamma for RF calculation
gam = 267.5221; %< rad /s /uT
%%% lineshape - using scalar delta
switch tissuepars.lineshape
    case 'SL'
        [G,w_loc] = SuperLorentzian_lineshape(T2s,delta,'interpzero');% seconds
    case 'Gaussian'
        [G,w_loc] = gauss_lineshape(T2s,delta);% seconds
end

%% Lambda matrix and C are time invariant
% Free pool transverse components - no off resonance here
La = [-R2f 0;0 -R2f];
% the rest
Lb = [-k*M0s-R1(1) k*M0f k*M0f 0;k*(1-f)*M0s -k*M0f-R1(2) 0 0;...
    k*f*M0s 0 -k*M0f-R1(3) 0;0 0 0 -R1(4)];
Lambda = blkdiag(La,Lb);
C = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

%% Get number of time points from b1_band
nt = size(b1_band, 1);

% Initialize evolution matrix
Xtilde_rf = eye(7);

%% Loop over all entries in b1_band vector
for tt = 1:nt
    % Get current b1 value
    b1_value = b1_band(tt);
    
    % Take the complex pulse components for free pool
    b1x = real(b1_value); % Real part (x-component) of complex pulse
    b1y = imag(b1_value); % Imaginary part (y-component) of complex pulse
    
    OmegaFree = gam*[0 0 -b1y; 0 0 b1x; b1y -b1x 0];

    
    % Semisolid pool
    switch nband
       case '1band'
           w1 = gam*abs(b1_value); 
       case '2band'
           w1 = gam*(max(b1_value)/sqrt(2)); 
    end

     w = pi*w1^2*G;
    
    % Calculate frequency offset effect
    if abs(delta) == 0
        D = 0; % w*D = 0 when on-resonance
    else
        D = 2*pi*delta/w_loc;
    end

    % OmegaSemi calculation depends on nband
    switch nband
        case '1band'
            OmegaSemi = [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]]; 
        case '2band'
            OmegaSemi = [[-2*w 0 0];[0 -2*w 0];[0 0 -2*w*D^2]]; 
    end
    
    Omega = blkdiag(OmegaFree, OmegaSemi);
    
    % Make overall evolution matrix
    Atilde = cat(1, [(Lambda+Omega) C], zeros(1,7));
    
    Xtilde_rf = expm(Atilde*dt) * Xtilde_rf;
end

% eigenvector methods
[V,D] = eig(Xtilde_rf);

% Find the eigenvalue with minimum absolute value
eigenvals = diag(D);
% Find index where eigenvalue = 1
[~, ind] = min(abs(eigenvals - 1));  % Find closest to 1

% Get the corresponding eigenvector
Mout = V(:, ind);

% [min_eval,ind] = min(abs(eigenvals));
% 
% % Verify it's actually zero (within machine precision)
% if min_eval > eps*100
%     error('No zero eigenvalue found. Minimum eigenvalue is %e. Check system setup.', min_eval);
% end


%% Calculate final magnetization
% Initialize the M_0 vector with 7 elements
% M0tilde = [0; 0; M0f; (1-f)*M0s; f*M0s; 0; 1];

% Calculate final state
% Mout = Xtilde_rf * M0tilde;

% Extract free water magnetization and normalize by M0f
Mz_out = abs(Mout(3));

% Debug output if requested
% if nargin > 4 && strcmp(varargin{1}, 'debug')
%     fprintf('Delta = %.1f Hz, nt = %d\n', delta, nt);
%     fprintf('Final Mz = %.6f, Mz_free = %.6f\n', abs(Mout(3)), Mz_free);
% end

end