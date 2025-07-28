% Z-spectrum simulation - Multiple B1 values for HC 1-band
%% pulse parameters
pulse_duration = 5; % saturation pulse duration
npoints = 100000; % samples in the shape
dt = pulse_duration/npoints; %dwell time

%% frequency offsets to sweep (positive range only)
offset_vec = linspace(100, 6e3, 31); % Hz (100 Hz to 6 kHz)

%% tissue parameters
tissuepars = init_tissue('hc');
tissuepars.lineshape = 'SL';

%% B1 values to test
B1_values = [3.45, 6.9, 6.9*sqrt(2), 6.9*2, 6.9*2*sqrt(2), 0]; % µT
B1_labels = {'3.45 µT', '6.9 µT', '6.9√2 µT', '13.8 µT', '13.8√2 µT', '0 µT (control)'};
colors = {'b-', 'r-', 'g-', 'm-', 'c-', 'k--'}; % Different colors/styles

%% Pre-allocate results
num_B1 = length(B1_values);
Mz_results = zeros(num_B1, length(offset_vec));

%% Simulation parameters (fixed for all B1 values)
nband = '1band';
shape = 'square';

%% Loop through each B1 value
for b = 1:num_B1
    B1_max = B1_values(b);
    
    for k = 1:numel(offset_vec)
        delta = offset_vec(k); % positive offset
        
        if B1_max == 0
            % Control case: no RF pulse - calculate equilibrium magnetization
            M0s = tissuepars.semi.M0;
            M0f = 1 - M0s; % Free water fraction  
            Mz_results(b, k) = M0f; % Equilibrium free water magnetization
        else
            % build the unit-max shape for THIS offset
            pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
            
            % scale to physical B1
            b1_single = B1_max * pulse_shape(:); % N×1, µT
            
            % distribute into the three-band table
            b1pulse = zeros(npoints, 3); % N×3
            b1pulse(:,3) = b1_single; % +Δ column only
            
            % set the Delta_Hz vector (always ±delta, 0)
            Delta_Hz = [-delta 0 +delta];
            
            % propagate once through Bloch-McConnell
            Mz_results(b, k) = ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
        end
    end
end

%% Normalize all results to M0 = 1
% The control case (0 µT) serves as M0 reference
M0_reference = Mz_results(end, :); % Last row is 0 µT case

% Normalize all cases including the control
for b = 1:num_B1
    Mz_results(b, :) = Mz_results(b, :) ./ M0_reference;
end

%% Plot all B1 cases on the same diagram
figure('DockControls','on');
hold on;

for b = 1:num_B1
    plot(offset_vec/1e3, Mz_results(b, :), colors{b}, 'LineWidth', 1.6, ...
         'DisplayName', B1_labels{b});
end

% Formatting
xlabel('\Delta (kHz)');
ylabel('$|M_{z,\mathrm{free}}|$ after RF pulse', 'Interpreter','latex');
title('HC 1-band Z-spectrum: B1 Power Comparison (9.4 T)');
legend('show', 'Location', 'best');
grid on;
hold off;

%% Save results

save('HC_1band_B1_comparison.mat', 'offset_vec', 'B1_values', 'Mz_results', ...
     'saturation_efficiency', 'tissuepars');