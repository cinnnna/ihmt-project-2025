%% Simultaneous Fitting of Hair Conditioner Datasets - CLEAN VERSION
% Fit ONLY the 7 hair conditioner datasets you specified

clear all; close all; clc;

%% Define EXACTLY the datasets you want (NO MORE, NO LESS)
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
    %{'Experimental Data/10Kegg_zspectra.mat', 'egg28dbdual', 'offset', '2band', 3.45*sqrt(2), '10K Egg 28dB dual'};
    
    % File 3: 0716eggwhite2.mat
    %{'Experimental Data/0716eggwhite2.mat', 'egg19dbpos', 'offset', '1band', 4.879, '0716 Egg 19dB pos'};
    %{'Experimental Data/0716eggwhite2.mat', 'egg19dbneg', 'offset', '1band', 4.879, '0716 Egg 19dB neg'};
    %{'Experimental Data/0716eggwhite2.mat', 'egg16dbdualart', 'offsetart', '2band', 4.879*sqrt(2), '0716 Egg 16dB dual art'};
    
    % File 4: 0716eggwhite.mat
    %{'Experimental Data/0716eggwhite.mat', 'egg25dbpos', 'offset', '1band', 6.9, '0716 Egg 25dB pos'};
    %{'Experimental Data/0716eggwhite.mat', 'egg25dbneg', 'offset', '1band', 6.9, '0716 Egg 25dB neg'};
    %{'Experimental Data/0716eggwhite.mat', 'egg22dbdual', 'offset', '2band', 9.75807, '0716 Egg 22dB dual'};
};


num_datasets_defined = length(datasets);
fprintf('=== DEFINED %d DATASETS FOR FITTING ===\n', num_datasets_defined );
for i = 1:length(datasets)
    fprintf('%d: %s\n', i, datasets{i}{6});
end
fprintf('==========================================\n');

%% Setup parameters
pulse_duration = 5;
npoints = 100000;
dt = pulse_duration/npoints;
shape = 'square';

%% Parameter bounds for hair conditioner
lb = [0.2, 1, 0.01, 0.1, 1, 0.1, 5e-6, 1];      % Lower bounds
ub = [2, 50, 0.3, 20, 100, 1, 50e-6, 50];       % Upper bounds
param_names = {'R1_free', 'R2_free', 'M0_semi', 'R1_semi', 'R1D_semi', 'f_semi', 'T2_semi', 'k'};

%% Initialize arrays - EXACTLY 7 elements
num_datasets = num_datasets_defined;
exp_data_array = cell(num_datasets, 1);
offset_data_array = cell(num_datasets, 1);
nband_array = cell(num_datasets, 1);
B1_max_array = cell(num_datasets, 1);
name_array = cell(num_datasets, 1);

%% Load each dataset individually
fprintf('\n=== LOADING PROCESS ===\n');
for idx = 1:num_datasets
    file_path = fullfile('/Users/heisenberg/Documents/MATLAB', datasets{idx}{1});
    data_field = datasets{idx}{2};
    offset_field = datasets{idx}{3};
    nband_val = datasets{idx}{4};
    B1_max_val = datasets{idx}{5};
    name_val = datasets{idx}{6};
    
    fprintf('Loading %d/%d: %s\n', idx, num_datasets, name_val);
    fprintf('  File: %s\n', datasets{idx}{1});
    fprintf('  Fields: %s, %s\n', data_field, offset_field);
    
    % Load file
    loaded_data = load(file_path);
    data_struct = loaded_data.data;
    
    % Extract specific fields
    if ~isfield(data_struct, data_field)
        available_fields = fieldnames(data_struct);
        fprintf('  Available fields: %s\n', strjoin(available_fields, ', '));
        error('Field "%s" not found in %s', data_field, file_path);
    end
    
    raw_exp_data = data_struct.(data_field);
    raw_offset_data = data_struct.(offset_field);
    
    % Filter frequency range
    if min(raw_offset_data) >= 0
        freq_mask = raw_offset_data >= 0 & raw_offset_data <= 10000;
    else
        freq_mask = raw_offset_data >= -6000 & raw_offset_data <= 6000;
    end
    
    % Store in arrays
    exp_data_array{idx} = raw_exp_data(freq_mask);
    offset_data_array{idx} = raw_offset_data(freq_mask);
    nband_array{idx} = nband_val;
    B1_max_array{idx} = B1_max_val;
    name_array{idx} = name_val;
    
    fprintf('  → Loaded %d points, range %.1f to %.1f kHz\n', ...
        sum(freq_mask), min(offset_data_array{idx})/1e3, max(offset_data_array{idx})/1e3);
end

%% Verification - MUST be exactly 7
fprintf('\n=== VERIFICATION ===\n');
fprintf('Loaded datasets: %d\n', length(exp_data_array));
fprintf('Expected: %d\n', num_datasets_defined);

if length(exp_data_array) ~= num_datasets_defined
    error('CRITICAL ERROR: Expected %d datasets but have %d!', num_datasets_defined, length(exp_data_array));
end

fprintf('SUCCESS: All %d datasets loaded correctly\n', num_datasets_defined);
for i = 1:num_datasets_defined
    fprintf('  %d: %s (%d points)\n', i, name_array{i}, length(exp_data_array{i}));
end

%% Initialize parameters for hair conditioner
tissuepars_init = init_tissue('hc');
tissuepars_init.lineshape = 'SL';

X_init = [tissuepars_init.free.R1, tissuepars_init.free.R2, tissuepars_init.semi.M0, ...
          tissuepars_init.semi.R1, tissuepars_init.semi.R1D, tissuepars_init.semi.f, ...
          tissuepars_init.semi.T2, tissuepars_init.k];

%% Test simulation with initial parameters
fprintf('\n=== TESTING INITIAL PARAMETERS ===\n');
try
    test_cost = calculate_total_cost(X_init, exp_data_array, offset_data_array, ...
        nband_array, B1_max_array, pulse_duration, npoints, dt, shape);
    fprintf('Initial cost: %.6f\n', test_cost);
catch test_err
    fprintf('ERROR in initial test: %s\n', test_err.message);
    return;
end

%% Optimization setup
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 300, ...
    'MaxFunctionEvaluations', 1500, 'FunctionTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, 'UseParallel', false);

cost_function = @(x) calculate_total_cost(x, exp_data_array, offset_data_array, ...
    nband_array, B1_max_array, pulse_duration, npoints, dt, shape);

%% Run optimization
fprintf('\n=== STARTING OPTIMIZATION ===\n');
fprintf('Fitting %d hair conditioner datasets simultaneously\n', num_datasets_defined);

tic;
[X_optimal, final_cost, exit_flag, opt_output] = fmincon(cost_function, X_init, ...
    [], [], [], [], lb, ub, [], options);
opt_time = toc;

%% Display results
fprintf('\n=== OPTIMIZATION RESULTS ===\n');
fprintf('Time: %.1f seconds\n', opt_time);
fprintf('Exit flag: %d\n', exit_flag);
fprintf('Final cost: %.6f\n', final_cost);
fprintf('Iterations: %d\n', opt_output.iterations);

fprintf('\nOPTIMAL PARAMETERS:\n');
for i = 1:length(param_names)
    fprintf('%-12s: %.6f\n', param_names{i}, X_optimal(i));
end

%% Calculate fit quality for each dataset
fprintf('\n=== INDIVIDUAL DATASET FITS ===\n');
figure('Name', 'Hair Conditioner Dataset Fits', 'Position', [100, 100, 1000, 700]);

total_points = 0;
total_sse = 0;
total_sst = 0;

for i = 1:num_datasets_defined
    % Simulate this specific dataset
    sim_result = simulate_single_dataset(X_optimal, offset_data_array{i}, ...
        pulse_duration, npoints, dt, nband_array{i}, shape, B1_max_array{i});
    
    % Calculate statistics
    errors = exp_data_array{i} - sim_result;
    sse = sum(errors.^2);
    rmse = sqrt(sse / length(errors));
    sst = sum((exp_data_array{i} - mean(exp_data_array{i})).^2);
    r_squared = 1 - sse/sst;
    
    total_points = total_points + length(errors);
    total_sse = total_sse + sse;
    total_sst = total_sst + sst;
    
    fprintf('%-20s: RMSE=%.4f, R²=%.3f\n', name_array{i}, rmse, r_squared);
    
    % Plot
    subplot(3, 3, i);
    plot(offset_data_array{i}/1e3, exp_data_array{i}, 'ro', 'MarkerSize', 3);
    hold on;
    plot(offset_data_array{i}/1e3, sim_result, 'b-', 'LineWidth', 1.5);
    xlabel('Δ (kHz)');
    ylabel('M_z');
    title(sprintf('%s (R²=%.3f)', name_array{i}, r_squared), 'FontSize', 8);
    grid on;
    hold off;
end

overall_rmse = sqrt(total_sse / total_points);
overall_r2 = 1 - total_sse/total_sst;

fprintf('\nOVERALL STATISTICS:\n');
fprintf('Total points: %d\n', total_points);
fprintf('Overall RMSE: %.6f\n', overall_rmse);
fprintf('Overall R²: %.4f\n', overall_r2);

%% Calculate and plot parameter correlation matrix
fprintf('\n=== Calculating Parameter Correlation Matrix ===\n');

try
    % Create combined prediction function for lsqcurvefit
    combined_predict_function = @(x, dummy) predict_all_datasets_combined(x, exp_data_array, offset_data_array, ...
        nband_array, B1_max_array, pulse_duration, npoints, dt, shape);
    
    % Combine all experimental data
    combined_exp_data = [];
    for i = 1:num_datasets_defined
        combined_exp_data = [combined_exp_data; exp_data_array{i}(:)];
    end
    
    % Use lsqcurvefit to get Jacobian
    [p_est, resnorm, residual, ~, ~, ~, J] = lsqcurvefit(combined_predict_function, X_optimal, ...
        1, combined_exp_data, [], [], optimset('Display','off','MaxFunEvals',1));
    
    % Calculate correlation matrix
    gamma = resnorm/(length(combined_exp_data) - length(p_est));
    pCov = gamma*inv(J'*J);
    pCorr = abs(pCov./sqrt(diag(pCov)*diag(pCov)'));
    pCorr = full(pCorr);
    
    % Plot correlation matrix
    figure('Name', 'Parameter Correlation Matrix');
    imagesc(pCorr);
    colorbar;
    title('Parameter Correlation Matrix');
    xlabel('Parameter Index');
    ylabel('Parameter Index');
    
    set(gca, 'XTick', 1:length(param_names), 'XTickLabel', param_names, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:length(param_names), 'YTickLabel', param_names);
    set(gca, 'XDir', 'reverse'); % Flip the x-axis as you wanted
    
    % Add correlation values as text
    for i = 1:length(param_names)
        for j = 1:length(param_names)
            if i ~= j
                text(j, i, sprintf('%.2f', pCorr(i,j)), 'HorizontalAlignment', 'center', ...
                    'Color', 'white', 'FontWeight', 'bold');
            end
        end
    end
    
    fprintf('Correlation matrix plotted successfully\n');
    
catch correlation_err
    fprintf('Warning: Correlation matrix calculation failed: %s\n', correlation_err.message);
end

%% Helper function for correlation calculation
function combined_sim_data = predict_all_datasets_combined(params, exp_data_array, offset_data_array, ...
    nband_array, B1_max_array, pulse_duration, npoints, dt, shape)
    
    combined_sim_data = [];
    for i = 1:length(exp_data_array)
        sim_data = simulate_single_dataset(params, offset_data_array{i}, ...
            pulse_duration, npoints, dt, nband_array{i}, shape, B1_max_array{i});
        combined_sim_data = [combined_sim_data; sim_data(:)];
    end
end


%% Save results
save('hair_conditioner_fitting_results.mat', 'datasets', 'X_optimal', 'final_cost', ...
    'exit_flag', 'opt_output', 'overall_rmse', 'overall_r2', 'param_names', 'name_array');

fprintf('\nResults saved to hair_conditioner_fitting_results.mat\n');
fprintf('Hair conditioner fitting completed successfully!\n');

%% SUPPORT FUNCTIONS

function total_cost = calculate_total_cost(params, exp_data_array, offset_data_array, ...
    nband_array, B1_max_array, pulse_duration, npoints, dt, shape)
    
    total_cost = 0;
    for i = 1:length(exp_data_array)
        sim_data = simulate_single_dataset(params, offset_data_array{i}, ...
            pulse_duration, npoints, dt, nband_array{i}, shape, B1_max_array{i});
        total_cost = total_cost + sum((exp_data_array{i} - sim_data).^2);
    end
end

function result = simulate_single_dataset(params, offset_vec, pulse_duration, npoints, dt, nband, shape, B1_max)
    
    % Convert parameters to tissue structure
    tissuepars.free.R1 = params(1);
    tissuepars.free.R2 = params(2);
    tissuepars.semi.M0 = params(3);
    tissuepars.semi.R1 = params(4);
    tissuepars.semi.R1D = params(5);
    tissuepars.semi.f = params(6);
    tissuepars.semi.T2 = params(7);
    tissuepars.k = params(8);
    tissuepars.lineshape = 'SL';
    
    result = zeros(size(offset_vec));
    
    for k = 1:length(offset_vec)
        current_offset = offset_vec(k);
        delta = abs(current_offset);
        sign_offset = sign(current_offset);
        
        if abs(current_offset) < eps
            sign_offset = 0;
        end
        
        % Generate pulse
        pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
        b1_amplitude = B1_max * pulse_shape(:);
        
        % Create 2-band pulse array
        b1pulse = zeros(npoints, 2);
        
        if delta == 0
            b1pulse(:,1) = b1_amplitude;
            b1pulse(:,2) = 0;
            Delta_Hz = [0, 0];
        else
            switch nband
                case '1band'
                    if sign_offset >= 0
                        b1pulse(:,1) = 0;
                        b1pulse(:,2) = b1_amplitude;
                        Delta_Hz = [0, delta];
                    else
                        b1pulse(:,1) = b1_amplitude;
                        b1pulse(:,2) = 0;
                        Delta_Hz = [-delta, 0];
                    end
                case '2band'
                    half_amplitude = b1_amplitude / sqrt(2);
                    b1pulse(:,1) = half_amplitude;
                    b1pulse(:,2) = half_amplitude;
                    Delta_Hz = [-delta, +delta];
            end
        end
        
        % Call your integration function
        result(k) = newDualcase_ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
    end
end