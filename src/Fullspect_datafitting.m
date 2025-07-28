%% MT Parameter Fitting using fmincon - MODIFIED FOR 2-BAND SYSTEM
% Fit simulation results to experimental data across full [-6000, 6000] Hz range

%% Load experimental data
%load('/Users/heisenberg/Documents/MATLAB/Experimental Data/6Khairconditioner_zspectra.mat');
load("/Users/heisenberg/Documents/MATLAB/Experimental Data/6Keggwhite_zspectra.mat");
%load("/Users/heisenberg/Documents/MATLAB/10Kegg_zspectra.mat");
%load("/Users/heisenberg/Documents/MATLAB/10KHC_zpectra.mat");

%% Setup simulation parameters (fixed)
pulse_duration = 5; % saturation pulse duration
npoints = 100000; % samples in the shape
dt = pulse_duration/npoints; %dwell time
nband = '1band'; % or '2band' depending on your experiment
shape = 'square';
B1_max = 6.9; % µT (peak B1) - adjust based on db25 power level

%% Initialize tissue parameters as vector X_in
% Get initial tissue parameters
tissuepars_init = init_tissue('random');
tissuepars_init.lineshape = 'SL';

% Convert tissue parameters to vector X_in
% Order: [R1_free, R2_free, M0_semi, R1_semi, R1D_semi, f_semi, T2_semi, k]
X_in = [tissuepars_init.free.R1, ...        % 1
        tissuepars_init.free.R2, ...        % 2
        tissuepars_init.semi.M0, ...        % 3
        tissuepars_init.semi.R1, ...        % 4
        tissuepars_init.semi.R1D, ...       % 5
        tissuepars_init.semi.f, ...         % 6
        tissuepars_init.semi.T2, ...        % 7
        tissuepars_init.k];                 % 8

%% Define parameter bounds for fmincon
% Lower bounds (reasonable physical limits) for eggwhite
lb = [0.2,    % R1_free (s^-1)
      1,      % R2_free (s^-1)  
      0.01,   % M0_semi (fraction)
      0.1,    % R1_semi (s^-1)
      1,      % R1D_semi (s^-1)
      0.1,    % f_semi (fraction) - FIXED: was 10, should be < 1
      5e-6,   % T2_semi (s)
      1];     % k (s^-1)

%Upper bounds for eggwhite
ub = [2,     % R1_free (s^-1)
      50,    % R2_free (s^-1)
      0.3,   % M0_semi (fraction)
      20,    % R1_semi (s^-1)
      100,   % R1D_semi (s^-1)
      1,     % f_semi (fraction)
      50e-6, % T2_semi (s)
      20];   % k (s^-1)

%% Setup experimental data for fitting - FULL RANGE
exp_data = data.eggdb25; % Using db25 experimental data
actual_offset_Hz = data.offset; % Data is already in Hz, no conversion needed

%% Extract experimental frequencies within FULL range (-6000 Hz to 6000 Hz)
range_indices = actual_offset_Hz >= -6000 & actual_offset_Hz <= 6000;
offset_vec = actual_offset_Hz(range_indices);
exp_data_interp = exp_data(range_indices);

% Display info
fprintf('Data used for fitting:\n');
fprintf('Frequency range: %.1f to %.1f kHz (%d points)\n', ...
    min(offset_vec)/1e3, max(offset_vec)/1e3, length(offset_vec));

% Test simulation with initial parameters
fprintf('Testing simulation with initial parameters...\n');
Mz_test = simulate_mt_spectrum_2band(X_in, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);
fprintf('Simulation range: %.6f to %.6f\n', min(Mz_test), max(Mz_test));
fprintf('First 5 simulation values: ');
for i = 1:min(5, length(Mz_test))
    fprintf('%.6f ', Mz_test(i));
end
fprintf('\n');

% Check initial error
initial_error = sum((exp_data_interp - Mz_test).^2);
fprintf('Initial error: %.10f\n', initial_error);
fprintf('========================\n\n');

%% Display frequency range info
fprintf('Using experimental frequencies within FULL -6000 to +6000 Hz range:\n');
fprintf('Number of frequencies: %d\n', length(offset_vec));
fprintf('Frequency range: %.1f to %.1f Hz\n', min(offset_vec), max(offset_vec));
fprintf('First few frequencies: ');
for i = 1:min(5, length(offset_vec))
    fprintf('%.1f ', offset_vec(i));
end
fprintf('Hz\n');

%% Setup optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ... 
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-6, ... 
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'UseParallel', true);

%% Run optimization
fprintf('Starting parameter optimization...\n');
fprintf('Initial parameters: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.2e, %.4f]\n', X_in);

tic;
[X_opt, fval, exitflag, output] = fmincon(@(x) errfun_2band(x, exp_data_interp, offset_vec, ...
    pulse_duration, npoints, dt, nband, shape, B1_max), ...
    X_in, [], [], [], [], lb, ub, [], options);
toc;

%% Display results
fprintf('\nOptimization completed!\n');
fprintf('Exit flag: %d\n', exitflag);
fprintf('Final error: %.6f\n', fval);
fprintf('Iterations: %d\n', output.iterations);

fprintf('\nOptimal parameters:\n');
fprintf('R1_free:  %.4f s^-1\n', X_opt(1));
fprintf('R2_free:  %.4f s^-1\n', X_opt(2));
fprintf('M0_semi:  %.4f\n', X_opt(3));
fprintf('R1_semi:  %.4f s^-1\n', X_opt(4));
fprintf('R1D_semi: %.4f s^-1\n', X_opt(5));
fprintf('f_semi:   %.4f\n', X_opt(6));
fprintf('T2_semi:  %.2e s\n', X_opt(7));
fprintf('k:        %.4f s^-1\n', X_opt(8));

%% Generate curves for comparison
Mz_initial = simulate_mt_spectrum_2band(X_in, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);
Mz_fitted = simulate_mt_spectrum_2band(X_opt, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);

figure('DockControls','on');
hold on;
plot(offset_vec/1e3, exp_data_interp, 'ro', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Experimental');
plot(offset_vec/1e3, Mz_initial, 'g--', 'LineWidth', 2, 'DisplayName', 'Initial simulation');
plot(offset_vec/1e3, Mz_fitted, 'b-', 'LineWidth', 2, 'DisplayName', 'Optimized simulation');
xlabel('\Delta (kHz)');
ylabel('Normalized |Mz_{free}|');
title('MT Parameter Fitting Results - Full Range');
legend('show', 'Location', 'best');
grid on;
xlim([-6.5, 6.5]); % Show full range
hold off;

%% Save results
tissuepars_fitted = vector_to_tissuepars(X_opt);
save('fitted_parameters_full_range.mat', 'X_opt', 'tissuepars_fitted', 'fval', 'output');

%% Error function for optimization - MODIFIED FOR 2-BAND SYSTEM
function err = errfun_2band(x, exp_data, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max)
    % Calculate simulated Z-spectrum using 2-band system
    Mz_sim = simulate_mt_spectrum_2band(x, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);
    
    % Calculate sum of squared errors
    err = sum((exp_data - Mz_sim).^2);
end

%% Function to simulate MT spectrum from parameter vector - MODIFIED FOR 2-BAND SYSTEM
function Mz_vec_norm = simulate_mt_spectrum_2band(x, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max)
    % Convert parameter vector to tissue parameters structure
    tissuepars = vector_to_tissuepars(x);
    
    % Pre-allocate results
    Mz_vec_norm = zeros(size(offset_vec));
    
    % Simulate each frequency offset using 2-band system
    for k = 1:length(offset_vec)
        current_offset = offset_vec(k);
        delta = abs(current_offset); % positive magnitude
        signOff = sign(current_offset); % +1 or -1
        
        % Handle zero offset case
        if abs(current_offset) < eps
            signOff = 0;
        end
        
        % Generate pulse
        pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
        b1_single = B1_max * pulse_shape(:);
        
        % Create 2-band b1pulse array
        b1pulse = zeros(npoints, 2); % N×2 for 2-band system
        
        % Apply pulse assignment logic for 2-band system
        if delta == 0
            % On-resonance case
            b1pulse(:,1) = b1_single; % 0 Hz band
            b1pulse(:,2) = 0;
            Delta_Hz = [0, 0];
        else
            switch nband
                case '1band'
                    if signOff >= 0
                        % Positive offset: [0, +delta] with pulse at +delta
                        b1pulse(:,1) = 0;
                        b1pulse(:,2) = b1_single;
                        Delta_Hz = [0, delta];
                    else
                        % Negative offset: [-delta, 0] with pulse at -delta  
                        b1pulse(:,1) = b1_single;
                        b1pulse(:,2) = 0;
                        Delta_Hz = [-delta, 0];
                    end
                case '2band'
                    % For 2band case: symmetrical ±Δ
                    half = b1_single / sqrt(2);
                    b1pulse(:,1) = half; % –Δ
                    b1pulse(:,2) = half; % +Δ
                    Delta_Hz = [-delta, +delta];
            end
        end
        
        % Run simulation using 2-band integrate function
        Mz_vec_norm(k) = Dualcase_ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
    end
end

%% Function to convert parameter vector to tissue parameters structure
function tissuepars = vector_to_tissuepars(x)
    tissuepars.free.R1 = x(1);
    tissuepars.free.R2 = x(2);
    tissuepars.semi.M0 = x(3);
    tissuepars.semi.R1 = x(4);
    tissuepars.semi.R1D = x(5);
    tissuepars.semi.f = x(6);
    tissuepars.semi.T2 = x(7);
    tissuepars.k = x(8);
    tissuepars.lineshape = 'SL'; % Changed to match your setup
end