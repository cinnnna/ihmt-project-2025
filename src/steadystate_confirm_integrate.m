%%% [Mz_free, Mz_evolution] = steadystate_confirm_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)
%
% Modified ssSPGR_ihMT_integrate to track steady-state convergence
%
% INPUTS:
% b1pulse = RF pulse, Mx3 array (M=#timepoints, 3=frequency
% bands). Units are uT
% dt = dwell time, sec
% Delta_Hz = 1x3 vector of frequencies of each band. Typically
% [-delta 0 delta]. Units Hz
% tissuepars = structure containing all tissue parameters. See
% init_tissue()
% varargin: optional 'track_convergence' flag to enable convergence tracking
%
% OUTPUTS:
% Mz_free = Longitudinal magnetization of free water (final value)
% Mz_evolution = Array of Mz_free values at each time step (if tracking enabled)
%
% (c) Shaihan Malik 2019. King's College London
% Modified to track steady-state convergence

function [Mz_free, Mz_evolution] = steadystate_confirm_integrate(b1pulse,dt,Delta_Hz,tissuepars,varargin)

% Check if convergence tracking is requested
track_convergence = false;
if nargin > 4 && any(strcmp(varargin, 'track_convergence'))
    track_convergence = true;
end

%%% Unpack tissue parameters
M0s = tissuepars.semi.M0; 
f = tissuepars.semi.f;
M0f = 1-M0s; 
R1 = [tissuepars.free.R1 tissuepars.semi.R1 tissuepars.semi.R1 tissuepars.semi.R1D];
R2f = tissuepars.free.R2;
T2s = tissuepars.semi.T2;
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
nt = size(b1pulse,1); %gives number of rows in b1pulse and shows how many discrete time-steps the pulse occupies
tau = dt*nt;

% Initialize the M_0 vector with 7 elements
M0tilde = [ 0 ; 0 ; M0f ; (1-f)*M0s ; f*M0s ; 0 ; 1 ];

%% Initialize tracking arrays if requested
if track_convergence
    Mz_evolution = zeros(nt, 1);
    Mout_current = M0tilde; % Start with initial state
else
    Mz_evolution = [];
end

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
    
    % apply matrix exponential and multiply into existing matrix
    Xtilde_rf = expm(Atilde*dt)*Xtilde_rf;
    
    % Track convergence if requested
    if track_convergence
        % Apply current transformation to get state at this time point
        Mout_current = Xtilde_rf * M0tilde;
        % Store the free water longitudinal magnetization
        Mz_evolution(tt) = abs(Mout_current(3));
    end
end

%% Now compile these and simulate longitudinal magnetization of free water
Mout = Xtilde_rf * M0tilde; % state right after the RF pulse
Mz_free = abs(Mout(3)); % the third element is the free-water longitudinal component

% Display convergence information if tracking was enabled
if track_convergence
    fprintf('Steady-state convergence analysis:\n');
    fprintf('Initial Mz_free: %.6f\n', abs(M0tilde(3)));
    fprintf('Final Mz_free: %.6f\n', Mz_free);
    fprintf('Change over pulse duration: %.6f\n', Mz_free - abs(M0tilde(3)));
    
    % Check convergence in the last 10% of the pulse
    last_10_percent = max(1, round(0.9*nt)):nt;
    if length(last_10_percent) > 1 %at least two points for calculating std
        final_variation = std(Mz_evolution(last_10_percent)); %low std means steady-state
        fprintf('Standard deviation in final 10%% of pulse: %.8f\n', final_variation);
        
        if final_variation < 0.005
            fprintf('✓ System appears to have reached steady state\n');
        else
            fprintf('⚠ System may not be at steady state yet\n');
        end
    end
    
    % Plot convergence if requested
    if any(strcmp(varargin, 'plot_convergence'))
        figure;
        time_axis = (1:nt) * dt;
        plot(time_axis, Mz_evolution, 'b-', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Mz_{free}');
        title('Magnetization Evolution During Pulse');
        grid on;
        
        % Add horizontal line at final value
        hold on;
        plot([time_axis(1), time_axis(end)], [Mz_free, Mz_free], 'r--', 'LineWidth', 1);
        legend('Mz evolution', 'Final value', 'Location', 'best');
        hold off;
    end
end

end