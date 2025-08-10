%% Simultaneous Fitting of Hair Conditioner Datasets - CLEAN VERSION
% Fit ONLY the 7 hair conditioner datasets you specified

clear all; close all; clc;

%% Define EXACTLY the datasets you want (NO MORE, NO LESS)
datasets = {
    % File 1: 6Keggwhite_zspectra.mat
    %{'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb19', 'offset', '1band', 4.879, '6K Egg dB19'};
    %{'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb25', 'offset', '1band', 6.9, '6K Egg dB25'};
    %{'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb31', 'offset', '1band', 3.45, '6K Egg dB31'};
    %{'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb37', 'offset', '1band', 1.7125, '6K Egg dB37'};
    %{'Experimental Data/6Keggwhite_zspectra.mat', 'eggdb43', 'offset', '1band', 0.8625, '6K Egg dB43'};
    %{'Experimental Data/6Keggwhite_zspectra.mat', 'eggdboff', 'offset', '1band', 0, '6K Egg dBoff'};
    
    % File 2: 10Kegg_zspectra.mat
    %{'Experimental Data/10Kegg_zspectra.mat', 'egg31dbsingle', 'offset', '1band', 3.45, '10K Egg 31dB single'};
    %{'Experimental Data/10Kegg_zspectra.mat', 'egg28dbdual', 'offset', '2band', 3.45*sqrt(2), '10K Egg 28dB dual'};
    
    % File 3: 0716eggwhite2.mat
    {'Experimental Data/0716eggwhite2.mat', 'egg19dbpos', 'offset', '1band', 4.879, '0716 Egg 19dB pos'};
    {'Experimental Data/0716eggwhite2.mat', 'egg19dbneg', 'offset', '1band', 4.879, '0716 Egg 19dB neg'};
    %{'Experimental Data/0716eggwhite2.mat', 'egg16dbdualart', 'offsetart', '2band', 4.879*sqrt(2), '0716 Egg 16dB dual art'};
    
    % File 4: 0716eggwhite.mat
    {'Experimental Data/0716eggwhite.mat', 'egg25dbpos', 'offset', '1band', 6.9, '0716 Egg 25dB pos'};
    {'Experimental Data/0716eggwhite.mat', 'egg25dbneg', 'offset', '1band', 6.9, '0716 Egg 25dB neg'};
    %{'Experimental Data/0716eggwhite.mat', 'egg22dbdual', 'offset', '2band', 9.75807, '0716 Egg 22dB dual'};
};

%% Pulse parameters
pulse_duration = 5; % seconds
npoints = 100000;
dt = pulse_duration/npoints;
shape = 'square';

%% Parameter bounds
%         R1f   R2f   M0s   R1s   R1D    f    T2s    k
lb =     [0.2,   1,  0.01,  0.1,   1,  0,  5e-6,  1];
ub =     [2,    50,  0.3,   20,  100,  0,   50e-6, 50];
param_names = {'R1_free', 'R2_free', 'M0_semi', 'R1_semi', 'R1D_semi', 'f_semi', 'T2_semi', 'k'};

%% Load experimental data
num_datasets = length(datasets);
exp_data = cell(num_datasets, 1);
offset_Hz = cell(num_datasets, 1);
nband = cell(num_datasets, 1);
B1_max = zeros(num_datasets, 1);
names = cell(num_datasets, 1);

for i = 1:num_datasets
    % Load data
    data = load(fullfile('/Users/heisenberg/Documents/MATLAB', datasets{i}{1}));
    exp_data{i} = data.data.(datasets{i}{2});
    offset_Hz{i} = data.data.(datasets{i}{3});
    nband{i} = datasets{i}{4};
    B1_max(i) = datasets{i}{5};
    names{i} = datasets{i}{6};
    
    % Filter frequency range
    if min(offset_Hz{i}) >= 0
        mask = offset_Hz{i} <= 6000;
    else
        mask = abs(offset_Hz{i}) <= 6000;
    end
    exp_data{i} = exp_data{i}(mask);
    offset_Hz{i} = offset_Hz{i}(mask);
end

%% Initial parameters from tissue defaults
tissuepars = init_tissue('hc');
X_init = [tissuepars.free.R1, tissuepars.free.R2, tissuepars.semi.M0, ...
          tissuepars.semi.R1, tissuepars.semi.R1D, tissuepars.semi.f, ...
          tissuepars.semi.T2, tissuepars.k];

%% Define cost function
cost_fun = @(x) calculate_cost(x, exp_data, offset_Hz, nband, B1_max, ...
                               pulse_duration, npoints, dt, shape);

%% Optimize
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 30);
[X_opt, final_cost] = fmincon(cost_fun, X_init, [], [], [], [], lb, ub, [], options);

%% Display optimized parameters
fprintf('\n=== OPTIMIZED PARAMETERS ===\n');
for i = 1:length(param_names)
    fprintf('%-12s: %.6f\n', param_names{i}, X_opt(i));
end

%% Calculate correlation matrix
% Calculate Jacobian numerically
epsilon = 1e-6;
n_params = length(X_opt);
J = [];

% Collect residuals at optimal parameters
for i = 1:num_datasets
    sim_opt = simulate_dataset(X_opt, offset_Hz{i}, pulse_duration, npoints, ...
                               dt, nband{i}, shape, B1_max(i));
    residuals_opt = exp_data{i} - sim_opt;
    
    % Calculate Jacobian for this dataset
    J_dataset = zeros(length(residuals_opt), n_params);
    for j = 1:n_params
        X_perturb = X_opt;
        X_perturb(j) = X_opt(j) + epsilon;
        sim_perturb = simulate_dataset(X_perturb, offset_Hz{i}, pulse_duration, ...
                                       npoints, dt, nband{i}, shape, B1_max(i));
        J_dataset(:,j) = (sim_opt - sim_perturb) / epsilon;
    end
    J = [J; J_dataset];
end

% Calculate correlation matrix
gamma = final_cost / (length(J) - n_params);
pCov = gamma * inv(J'*J);
pCorr = abs(pCov ./ sqrt(diag(pCov) * diag(pCov)'));

% Plot correlation matrix
figure;
imagesc(pCorr);
colorbar;
title('Parameter Correlation Matrix');
set(gca, 'XTick', 1:n_params, 'XTickLabel', param_names);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names);
xtickangle(45);
set(gca, 'XDir', 'reverse');
saveas(gcf, 'correlation_matrix.fig');

%% Plot fits and calculate statistics
figure('Position', [100, 100, 1200, 800]);
total_sse = 0;
total_sst = 0;
total_points = 0;

for i = 1:num_datasets
    % Simulate with optimized parameters
    sim_data = simulate_dataset(X_opt, offset_Hz{i}, pulse_duration, npoints, ...
                               dt, nband{i}, shape, B1_max(i));
    
    % Calculate statistics
    residuals = exp_data{i} - sim_data;
    sse = sum(residuals.^2);
    sst = sum((exp_data{i} - mean(exp_data{i})).^2);
    rmse = sqrt(sse / length(residuals));
    nrmse = rmse / (max(exp_data{i}) - min(exp_data{i}));
    r2 = 1 - sse/sst;
    
    total_sse = total_sse + sse;
    total_sst = total_sst + sst;
    total_points = total_points + length(residuals);
    
    % Plot
    subplot(2, 3, i);
    plot(offset_Hz{i}/1e3, exp_data{i}, 'ro', 'MarkerSize', 4, 'DisplayName', 'Experimental');
    hold on;
    plot(offset_Hz{i}/1e3, sim_data, 'b-', 'LineWidth', 2, 'DisplayName', 'Fit');
    xlabel('Frequency offset (kHz)');
    ylabel('M_z/M_0');
    title(sprintf('%s\nRMSE=%.4f, NRMSE=%.2f%%, R²=%.3f', ...
                  names{i}, rmse, nrmse*100, r2));
    legend('Location', 'best');
    grid on;
    xlim([min(offset_Hz{i})/1e3, max(offset_Hz{i})/1e3]);
    ylim([0, 1.1]);
end

% Overall statistics
overall_rmse = sqrt(total_sse / total_points);
overall_r2 = 1 - total_sse/total_sst;
fprintf('\n=== OVERALL FIT QUALITY ===\n');
fprintf('Overall RMSE: %.6f\n', overall_rmse);
fprintf('Overall R²: %.4f\n', overall_r2);



%% === HELPER FUNCTIONS ===

function cost = calculate_cost(params, exp_data, offset_Hz, nband, B1_max, ...
                              pulse_duration, npoints, dt, shape)
    cost = 0;
    for i = 1:length(exp_data)
        sim = simulate_dataset(params, offset_Hz{i}, pulse_duration, npoints, ...
                              dt, nband{i}, shape, B1_max(i));
        cost = cost + sum((exp_data{i} - sim).^2);
    end
end

function sim_data = simulate_dataset(params, offset_vec, pulse_duration, ...
                                    npoints, dt, nband, shape, B1_max)
    % Create tissue parameters structure
    tissuepars.free.R1 = params(1);
    tissuepars.free.R2 = params(2);
    tissuepars.semi.M0 = params(3);
    tissuepars.semi.R1 = params(4);
    tissuepars.semi.R1D = params(5);
    tissuepars.semi.f = params(6);
    tissuepars.semi.T2 = params(7);
    tissuepars.k = params(8);
    tissuepars.lineshape = 'Gaussian';
    
    sim_data = zeros(size(offset_vec));
    
    for k = 1:length(offset_vec)
        delta = abs(offset_vec(k));
        
        % Generate pulse shape
        pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
        b1_band = B1_max * pulse_shape(:);
        
        % Simulate
        sim_data(k) = new_Dualcase_ssSPGR_ihMT_integrate(b1_band, dt, delta, tissuepars, nband);
    end
end
