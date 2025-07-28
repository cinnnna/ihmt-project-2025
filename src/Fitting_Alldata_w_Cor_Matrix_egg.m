%% Simultaneous Fitting of Specific Experimental Datasets - FIXED
% Fit multiple experimental datasets with different B1 and nband settings

clear; close all;

%% Define datasets with their specific parameters
datasets = {
    % File 1: 6Keggwhite_zspectra.mat
    {'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb19', 'offset', '1band', 4.879, '6K Egg dB19'};
    {'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb25', 'offset', '1band', 6.9, '6K Egg dB25'};
    {'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb31', 'offset', '1band', 3.45, '6K Egg dB31'};
    {'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb37', 'offset', '1band', 1.7125, '6K Egg dB37'};
    {'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb43', 'offset', '1band', 0.8625, '6K Egg dB43'};
    {'Experimental Data/6Keggwhite_zspectra.mat', 'eggdboff', 'offset', '1band', 0, '6K Egg dBoff'};
    
    % File 2: 10Kegg_zspectra.mat
    {'Experimental Data/10Kegg_zspectra.mat', 'egg31dbsingle', 'offset', '1band', 3.45, '10K Egg 31dB single'};
    {'Experimental Data/10Kegg_zspectra.mat', 'egg28dbdual', 'offset', '2band', 3.45*sqrt(2), '10K Egg 28dB dual'};
    
    % File 3: 0716eggwhite2.mat
    {'Experimental Data/0716eggwhite2.mat', 'egg19dbpos', 'offset', '1band', 4.879, '0716 Egg 19dB pos'};
    {'Experimental Data/0716eggwhite2.mat', 'egg19dbneg', 'offset', '1band', 4.879, '0716 Egg 19dB neg'};
    {'Experimental Data/0716eggwhite2.mat', 'egg16dbdualart', 'offsetart', '2band', 4.879*sqrt(2), '0716 Egg 16dB dual art'};
    
    % File 4: 0716eggwhite.mat
    {'Experimental Data/0716eggwhite.mat', 'egg25dbpos', 'offset', '1band', 6.9, '0716 Egg 25dB pos'};
    {'Experimental Data/0716eggwhite.mat', 'egg25dbneg', 'offset', '1band', 6.9, '0716 Egg 25dB neg'};
    {'Experimental Data/0716eggwhite.mat', 'egg22dbdual', 'offset', '2band', 9.75807, '0716 Egg 22dB dual'};
};

%% Setup simulation parameters (fixed)
pulse_duration = 5; % saturation pulse duration
npoints = 100000; % samples in the shape
dt = pulse_duration/npoints; %dwell time
shape = 'square';

%% Parameter bounds
% Lower bounds: [R1_free, R2_free, M0_semi, R1_semi, R1D_semi, f_semi, T2_semi, k]
lb = [0.2,    % R1_free (s^-1)
      1,      % R2_free (s^-1)  
      0.01,   % M0_semi (fraction)
      0.1,    % R1_semi (s^-1)
      1,      % R1D_semi (s^-1)
      0.1,    % f_semi (fraction)
      5e-6,   % T2_semi (s)
      1];     % k (s^-1)

% Upper bounds
ub = [2,     % R1_free (s^-1)
      50,    % R2_free (s^-1)
      0.3,   % M0_semi (fraction)
      20,    % R1_semi (s^-1)
      100,   % R1D_semi (s^-1)
      1,     % f_semi (fraction)
      50e-6, % T2_semi (s)
      50];   % k (s^-1)

%% Load and prepare all datasets
fprintf('=== Loading All Datasets ===\n');
num_datasets = length(datasets);
all_exp_data = {};
all_offset_data = {};
all_nband = {};
all_B1_max = {};
dataset_names = {};

for i = 1:num_datasets
    file_path = fullfile('/Users/heisenberg/Documents/MATLAB', datasets{i}{1});
    data_field = datasets{i}{2};
    offset_field = datasets{i}{3};
    nband = datasets{i}{4};
    B1_max = datasets{i}{5};
    name = datasets{i}{6};
    
    try
        fprintf('Loading dataset %d: %s\n', i, name);
        load(file_path);
        
        % Extract data using dynamic field names
        exp_data = data.(data_field);
        offset_data = data.(offset_field);
        
        % Use all available data within reasonable range (exclude very high frequencies)
        % Most of your data appears to be positive offsets from 0 to 100kHz
        % Exclude points above 78kHz to avoid spline interpolation issues
        if min(offset_data) >= 0
            % Positive-only offset data (0 to 78kHz to avoid spline issues)
            range_indices = offset_data >= 0 & offset_data <= 10000; % 0 to 78kHz
        else
            % Full range data (-6kHz to +6kHz)
            range_indices = offset_data >= -6000 & offset_data <= 6000;
        end
        
        all_exp_data{i} = exp_data(range_indices);
        all_offset_data{i} = offset_data(range_indices);
        all_nband{i} = nband;
        all_B1_max{i} = B1_max;
        dataset_names{i} = name;
        
        fprintf('  Loaded %d data points, range: %.1f to %.1f kHz\n', ...
            sum(range_indices), min(all_offset_data{i})/1e3, max(all_offset_data{i})/1e3);
        
    catch load_err
        error('Failed to load dataset %d (%s): %s', i, name, load_err.message);
    end
end

%% Initialize parameters
tissuepars_init = init_tissue('random');
tissuepars_init.lineshape = 'Gaussian';  % Changed to Gaussian for egg white

X_init = [tissuepars_init.free.R1, ...
          tissuepars_init.free.R2, ...
          tissuepars_init.semi.M0, ...
          tissuepars_init.semi.R1, ...
          tissuepars_init.semi.R1D, ...
          tissuepars_init.semi.f, ...
          tissuepars_init.semi.T2, ...
          tissuepars_init.k];

param_names = {'R1_free', 'R2_free', 'M0_semi', 'R1_semi', 'R1D_semi', 'f_semi', 'T2_semi', 'k'};

%% Setup optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ... 
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-8, ... 
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-8, ...
    'UseParallel', false);

%% Define combined cost function for simultaneous fitting
combined_cost_function = @(x) compute_combined_cost(x, all_exp_data, all_offset_data, ...
    all_nband, all_B1_max, pulse_duration, npoints, dt, shape);

%% Test initial parameters
fprintf('\n=== Testing Initial Parameters ===\n');
initial_cost = combined_cost_function(X_init);
fprintf('Initial combined cost: %.8f\n', initial_cost);

%% Run simultaneous optimization
fprintf('\n=== Starting Simultaneous Optimization ===\n');
fprintf('Fitting %d datasets simultaneously...\n', num_datasets);
fprintf('Initial parameters: ');
for i = 1:length(X_init)
    fprintf('%.4f ', X_init(i));
end
fprintf('\n');

tic;
[X_opt, fval, exitflag, output] = fmincon(combined_cost_function, X_init, ...
    [], [], [], [], lb, ub, [], options);
fitting_time = toc;

%% Display optimization results
fprintf('\n=== OPTIMIZATION RESULTS ===\n');
fprintf('Optimization completed in %.1f seconds\n', fitting_time);
fprintf('Exit flag: %d\n', exitflag);
fprintf('Final combined cost: %.8f\n', fval);
fprintf('Cost reduction: %.2f%%\n', 100*(initial_cost - fval)/initial_cost);
fprintf('Iterations: %d\n', output.iterations);
fprintf('Function evaluations: %d\n', output.funcCount);

%% Display best-fit parameters
fprintf('\nBEST-FIT TISSUE PARAMETERS:\n');
fprintf('Parameter\t\tOptimal Value\n');
fprintf('--------------------------------\n');
for i = 1:length(param_names)
    if i == 7 % T2_semi in scientific notation
        fprintf('%-12s\t\t%.2e s\n', param_names{i}, X_opt(i));
    else
        fprintf('%-12s\t\t%.6f\n', param_names{i}, X_opt(i));
    end
end

%% Calculate individual dataset fits and errors
fprintf('\nINDIVIDUAL DATASET FIT QUALITY:\n');
fprintf('Dataset\t\t\t\t\tRMSE\t\tR²\t\tMax Error\n');
fprintf('----------------------------------------------------------------\n');

total_points = 0;
total_sse = 0;
total_sst = 0;

figure('Name', 'Individual Dataset Fits', 'Position', [100, 100, 1200, 800]);

for i = 1:num_datasets
    % Simulate this dataset with optimal parameters
    sim_data = simulate_mt_spectrum_2band(X_opt, all_offset_data{i}, ...
        pulse_duration, npoints, dt, all_nband{i}, shape, all_B1_max{i});
    
    % Calculate fit statistics
    residuals = all_exp_data{i} - sim_data;
    sse = sum(residuals.^2);
    rmse = sqrt(sse / length(residuals));
    sst = sum((all_exp_data{i} - mean(all_exp_data{i})).^2);
    r_squared = 1 - sse/sst;
    max_error = max(abs(residuals));
    
    % Accumulate for overall statistics
    total_points = total_points + length(residuals);
    total_sse = total_sse + sse;
    total_sst = total_sst + sst;
    
    fprintf('%-25s\t\t%.6f\t%.4f\t\t%.6f\n', dataset_names{i}, rmse, r_squared, max_error);
    
    % Plot individual fits
    subplot(4, 4, i);
    plot(all_offset_data{i}/1e3, all_exp_data{i}, 'ro', 'MarkerSize', 3, 'DisplayName', 'Experimental');
    hold on;
    plot(all_offset_data{i}/1e3, sim_data, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Fitted');
    xlabel('\Delta (kHz)');
    ylabel('M_z');
    title(sprintf('%s (R²=%.3f)', dataset_names{i}, r_squared), 'FontSize', 8);
    legend('Location', 'best', 'FontSize', 6);
    grid on;
    hold off;
end

% Overall statistics
overall_rmse = sqrt(total_sse / total_points);
overall_r_squared = 1 - total_sse/total_sst;

fprintf('\nOVERALL FIT STATISTICS:\n');
fprintf('Total data points: %d\n', total_points);
fprintf('Overall RMSE: %.6f\n', overall_rmse);
fprintf('Overall R²: %.4f\n', overall_r_squared);

%% Calculate parameter uncertainties using combined Jacobian
fprintf('\n=== Calculating Parameter Uncertainties ===\n');

try
    % Create combined prediction function for lsqcurvefit
    combined_predict_function = @(x, dummy) predict_all_datasets(x, all_exp_data, all_offset_data, ...
        all_nband, all_B1_max, pulse_duration, npoints, dt, shape);
    
    % Combine all experimental data
    combined_exp_data = [];
    for i = 1:num_datasets
        combined_exp_data = [combined_exp_data; all_exp_data{i}(:)];
    end
    
    % Use lsqcurvefit to get Jacobian
    [p_est, resnorm, residual, ~, ~, ~, J] = lsqcurvefit(combined_predict_function, X_opt, ...
        1, combined_exp_data, [], [], optimset('Display','off','MaxFunEvals',1));
    
    % Calculate parameter uncertainties
    gamma = resnorm/(length(combined_exp_data) - length(p_est));
    pCov = gamma*inv(J'*J);
    pSD = sqrt(diag(pCov))';
    pCV = 100*pSD./abs(p_est);
    pCorr = abs(pCov./sqrt(diag(pCov)*diag(pCov)'));
    
    % Convert to full arrays
    p_est = full(p_est);
    pSD = full(pSD);
    pCV = full(pCV);
    pCorr = full(pCorr);
    
    % Display uncertainties
    fprintf('\nPARAMETER UNCERTAINTIES:\n');
    fprintf('Parameter\t\tValue ± Std Dev\t\tCV (%%)\n');
    fprintf('------------------------------------------------\n');
    for i = 1:length(param_names)
        fprintf('%-12s\t\t%.6f ± %.6f\t%.2f\n', param_names{i}, p_est(i), pSD(i), pCV(i));
    end
    
    % Plot correlation matrix
    figure('Name', 'Parameter Correlation Matrix');
    imagesc(pCorr);
    colorbar;
    title('Parameter Correlation Matrix');
    xlabel('Parameter Index');
    ylabel('Parameter Index');
    
    set(gca, 'XTick', 1:length(param_names), 'XTickLabel', param_names, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:length(param_names), 'YTickLabel', param_names);
    set(gca, 'XDir', 'reverse'); % Flip the x-axis
    
    % Add correlation values as text
    for i = 1:length(param_names)
        for j = 1:length(param_names)
            if i ~= j
                text(j, i, sprintf('%.2f', pCorr(i,j)), 'HorizontalAlignment', 'center', ...
                    'Color', 'white', 'FontWeight', 'bold');
            end
        end
    end
    
catch uncertainty_err
    fprintf('Warning: Parameter uncertainty calculation failed: %s\n', uncertainty_err.message);
    pSD = NaN(size(X_opt));
    pCV = NaN(size(X_opt));
    pCorr = NaN(length(X_opt), length(X_opt));
end

%% Save comprehensive results
save('simultaneous_fitting_results.mat', ...
    'datasets', 'dataset_names', 'param_names', 'X_opt', 'fval', 'exitflag', 'output', ...
    'pSD', 'pCV', 'pCorr', 'overall_rmse', 'overall_r_squared', 'fitting_time');

fprintf('\nResults saved to simultaneous_fitting_results.mat\n');
fprintf('Simultaneous fitting completed successfully!\n');

%% SUPPORTING FUNCTIONS

function total_cost = compute_combined_cost(x, all_exp_data, all_offset_data, all_nband, all_B1_max, ...
    pulse_duration, npoints, dt, shape)
    % Compute combined cost across all datasets
    total_cost = 0;
    
    for i = 1:length(all_exp_data)
        % Simulate this dataset
        sim_data = simulate_mt_spectrum_2band(x, all_offset_data{i}, ...
            pulse_duration, npoints, dt, all_nband{i}, shape, all_B1_max{i});
        
        % Add to total cost (sum of squared errors)
        total_cost = total_cost + sum((all_exp_data{i} - sim_data).^2);
    end
end

function combined_sim_data = predict_all_datasets(x, all_exp_data, all_offset_data, all_nband, all_B1_max, ...
    pulse_duration, npoints, dt, shape)
    % Predict all datasets for uncertainty calculation
    combined_sim_data = [];
    
    for i = 1:length(all_exp_data)
        sim_data = simulate_mt_spectrum_2band(x, all_offset_data{i}, ...
            pulse_duration, npoints, dt, all_nband{i}, shape, all_B1_max{i});
        combined_sim_data = [combined_sim_data; sim_data(:)];
    end
end

function Mz_vec_norm = simulate_mt_spectrum_2band(x, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max)
    % Simulate MT spectrum using 2-band system - WITHOUT problematic error handling
    tissuepars = vector_to_tissuepars(x);
    Mz_vec_norm = zeros(size(offset_vec));
    
    for k = 1:length(offset_vec)
        current_offset = offset_vec(k);
        delta = abs(current_offset);
        signOff = sign(current_offset);
        
        if abs(current_offset) < eps
            signOff = 0;
        end
        
        pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
        b1_single = B1_max * pulse_shape(:);
        
        b1pulse = zeros(npoints, 2);
        
        if delta == 0
            b1pulse(:,1) = b1_single;
            b1pulse(:,2) = 0;
            Delta_Hz = [0, 0];
        else
            switch nband
                case '1band'
                    if signOff >= 0
                        b1pulse(:,1) = 0;
                        b1pulse(:,2) = b1_single;
                        Delta_Hz = [0, delta];
                    else
                        b1pulse(:,1) = b1_single;
                        b1pulse(:,2) = 0;
                        Delta_Hz = [-delta, 0];
                    end
                case '2band'
                    half = b1_single / sqrt(2);
                    b1pulse(:,1) = half;
                    b1pulse(:,2) = half;
                    Delta_Hz = [-delta, +delta];
            end
        end
        
        % Direct integration without error handling - same as your working code
        Mz_vec_norm(k) = Dualcase_ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
    end
end

function tissuepars = vector_to_tissuepars(x)
    % Convert parameter vector to tissue parameters structure
    tissuepars.free.R1 = x(1);
    tissuepars.free.R2 = x(2);
    tissuepars.semi.M0 = x(3);
    tissuepars.semi.R1 = x(4);
    tissuepars.semi.R1D = x(5);
    tissuepars.semi.f = x(6);
    tissuepars.semi.T2 = x(7);
    tissuepars.k = x(8);
    tissuepars.lineshape = 'Gaussian';  % Changed to Gaussian for egg white
end