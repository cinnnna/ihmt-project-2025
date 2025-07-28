%% MT Parameter Fitting using fmincon
% Fit simulation results to experimental data
%% CONSISTENCY TEST: Generate synthetic "experimental" data
%% Setup simulation parameters (fixed)
pulse_duration = 5; % saturation pulse duration
npoints = 100000; % samples in the shape
dt = pulse_duration/npoints; %dwell time
nband = '1band'; % or '2band' depending on your experiment
shape = 'square';
B1_max = 6.9; % µT (peak B1) - adjust based on db25 power level

% Define known tissue parameters (ground truth)
tissuepars_true = init_tissue('hc');
tissuepars_true.lineshape = 'SL';

% Convert to parameter vector (ground truth)
X_true = [tissuepars_true.free.R1, tissuepars_true.free.R2, tissuepars_true.semi.M0, ...
          tissuepars_true.semi.R1, tissuepars_true.semi.R1D, tissuepars_true.semi.f, ...
          tissuepars_true.semi.T2, tissuepars_true.k];

fprintf('=== CONSISTENCY TEST ===\n');
fprintf('Ground truth parameters:\n');
fprintf('R1_free: %.4f, R2_free: %.4f, M0_semi: %.4f, R1_semi: %.4f\n', X_true(1:4));
fprintf('R1D_semi: %.4f, f_semi: %.4f, T2_semi: %.2e, k: %.4f\n', X_true(5:8));

% Define frequency range for synthetic data
offset_vec_synth = linspace(100, 6000, 30); % 30 points from 100 Hz to 6 kHz

% Generate synthetic Z-spectrum data
fprintf('Generating synthetic experimental data...\n');
synth_data = simulate_mt_spectrum(X_true, offset_vec_synth, pulse_duration, npoints, dt, nband, shape, B1_max);

fprintf('Synthetic data generated: %d points\n', length(synth_data));
fprintf('========================\n\n');

%% Load experimental data


%% Initialize tissue parameters as vector X_in
% Get initial tissue parameters
tissuepars_init = init_tissue('hc');
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
% Lower bounds (reasonable physical limits)
lb = [0.1,    % R1_free (s^-1)
      1,      % R2_free (s^-1)  
      0.01,   % M0_semi (fraction)
      0.1,    % R1_semi (s^-1)
      1,      % R1D_semi (s^-1)
      0,      % f_semi (fraction)
      1e-6,   % T2_semi (s)
      1];     % k (s^-1)

% Upper bounds
ub = [10,     % R1_free (s^-1)
      100,    % R2_free (s^-1)
      0.99,   % M0_semi (fraction)
      100,    % R1_semi (s^-1)
      1000,   % R1D_semi (s^-1)
      1,      % f_semi (fraction)
      1e-3,   % T2_semi (s)
      1000];  % k (s^-1)
%% Use synthetic data as "experimental" data for consistency test
exp_data = synth_data;
actual_offset_Hz = offset_vec_synth;

%% Extract experimental frequencies within desired range (100 Hz to 6000 Hz)
range_indices = actual_offset_Hz >= 100 & actual_offset_Hz <= 6000;
offset_vec = actual_offset_Hz(range_indices);
exp_data_interp = exp_data(range_indices);


% Display info
fprintf('Data used for fitting:\n');
fprintf('Frequency range: %.1f to %.1f kHz (%d points)\n', ...
    min(offset_vec)/1e3, max(offset_vec)/1e3, length(offset_vec));

% Test simulation with initial parameters
fprintf('Testing simulation with initial parameters...\n');
Mz_test = simulate_mt_spectrum(X_in, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);
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
fprintf('Using experimental frequencies within 100-6000 Hz range:\n');
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
    'MaxIterations', 5000, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'UseParallel', false);

%% Run optimization
fprintf('Starting parameter optimization...\n');
fprintf('Initial parameters: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.2e, %.4f]\n', X_in);

tic;
[X_opt, fval, exitflag, output] = fmincon(@(x) errfun(x, exp_data_interp, offset_vec, ...
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
Mz_initial = simulate_mt_spectrum(X_in, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);
Mz_fitted = simulate_mt_spectrum(X_opt, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);

figure('DockControls','on');
hold on;
plot(offset_vec/1e3, exp_data_interp, 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Experimental');
plot(offset_vec/1e3, Mz_initial, 'g--', 'LineWidth', 2, 'DisplayName', 'Initial simulation');
plot(offset_vec/1e3, Mz_fitted, 'b-', 'LineWidth', 2, 'DisplayName', 'Optimized simulation');
xlabel('\Delta (kHz)');
ylabel('Normalized |Mz_{free}|');
title('MT Parameter Fitting Results');
legend('show', 'Location', 'best');
grid on;
hold off;

%% Save results
tissuepars_fitted = vector_to_tissuepars(X_opt);
save('fitted_parameters.mat', 'X_opt', 'tissuepars_fitted', 'fval', 'output');

%% Error function for optimization
function err = errfun(x, exp_data, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max)
    % Calculate simulated Z-spectrum
    Mz_sim = simulate_mt_spectrum(x, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max);
    
    % Calculate sum of squared errors
    err = sum((exp_data - Mz_sim).^2);

end

%% Function to simulate MT spectrum from parameter vector
function Mz_vec_norm = simulate_mt_spectrum(x, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max)
    % Convert parameter vector to tissue parameters structure
    tissuepars = vector_to_tissuepars(x);
    
    % Pre-allocate results
    Mz_vec_norm = zeros(size(offset_vec));
    
    % Simulate each frequency offset using eigenvector method
    for k = 1:length(offset_vec)
        delta = offset_vec(k);
        
        % Avoid delta = 0 which causes spline interpolation issues
        if abs(delta) < 10  % If delta is effectively zero
            delta = sign(delta) * 10;  % Set to small non-zero value
            if delta == 0
                delta = 10;  % Handle exact zero case
            end
        end
        
        % Generate pulse
        pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
        b1_single = B1_max * pulse_shape(:);
        
        % Create b1pulse array
        b1pulse = zeros(npoints, 3);
        switch nband
            case '1band'
                b1pulse(:,3) = b1_single; % +Δ column only
            case '2band'
                half = b1_single / sqrt(2);
                b1pulse(:,1) = half; % -Δ
                b1pulse(:,3) = half; % +Δ
        end
        
        % Set frequency vector
        Delta_Hz = [-delta 0 +delta];
        
        % Run simulation using eigenvector method
        % (returns already normalized to M0f = 1)
        Mz_vec_norm(k) = simplified_ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
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
    tissuepars.lineshape = 'Gaussian';
end


%% CONSISTENCY CHECK: Compare recovered vs ground truth parameters
fprintf('\n=== CONSISTENCY CHECK RESULTS ===\n');
fprintf('Parameter comparison (Ground Truth vs Recovered):\n');
param_names = {'R1_free', 'R2_free', 'M0_semi', 'R1_semi', 'R1D_semi', 'f_semi', 'T2_semi', 'k'};

fprintf('%-12s %12s %12s %12s\n', 'Parameter', 'Ground Truth', 'Recovered', 'Rel. Error %');
fprintf('%-12s %12s %12s %12s\n', '----------', '-----------', '---------', '-----------');

for i = 1:8
    rel_error = 100 * abs(X_opt(i) - X_true(i)) / abs(X_true(i));
    fprintf('%-12s %12.6f %12.6f %12.2f\n', param_names{i}, X_true(i), X_opt(i), rel_error);
end

% Overall assessment
max_error = max(100 * abs(X_opt - X_true) ./ abs(X_true));
fprintf('\nMaximum relative error: %.2f%%\n', max_error);

if max_error < 1.0
    fprintf('✅ EXCELLENT: All parameters recovered within 1%%\n');
elseif max_error < 5.0
    fprintf('✅ GOOD: All parameters recovered within 5%%\n');
elseif max_error < 10.0
    fprintf('⚠️  ACCEPTABLE: Parameters recovered within 10%%\n');
else
    fprintf('❌ POOR: Some parameters have >10%% error\n');
end
fprintf('=====================================\n');