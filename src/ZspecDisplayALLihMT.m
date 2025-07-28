
% Z-spectrum simulation - Combined 1-band and 2-band comparison

%% pulse parameters
pulse_duration = 5; % saturation pulse duration
npoints = 100000; % samples in the shape
 
dt = pulse_duration/npoints; %dwell time

%% frequency offsets to sweep (positive range only)
offset_vec = linspace(100, 6e3, 31); % Hz (100 Hz to 5.5 kHz)

%% tissue parameters
tissuepars = init_tissue('hc');
tissuepars.lineshape = 'SL';

%% Pre-allocate results for both cases
Mz_vec_1band = zeros(size(offset_vec));
Mz_vec_2band = zeros(size(offset_vec));

%% Simulate 1-band case
%fprintf('Simulating 1-band case...\n');
nband = '1band';
shape = 'square';
B1_max = 6.9; % µT (peak B1)
for k = 1:numel(offset_vec)
    delta = offset_vec(k); % positive offset
    signOff = 1; % always positive for this range
    
    % build the unit-max shape for THIS offset
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
    
    % scale to physical B1
    b1_single = B1_max * pulse_shape(:); % N×1, µT
    
    % distribute into the three-band table
    b1pulse = zeros(npoints, 3); % N×3
    b1pulse(:,3) = b1_single; % +Δ column only (since signOff >= 0)
    
    % set the Delta_Hz vector (always ±delta, 0)
    Delta_Hz = [-delta 0 +delta];
    
    % propagate once through Bloch-McConnell
    Mz_vec_1band(k) = ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
    
end

%% Simulate 2-band case
%fprintf('Simulating 2-band case...\n');
nband = '2band';
shape = 'square';
B1_max = 6.9*sqrt(2) ;% µT (peak B1)
for k = 1:numel(offset_vec)
    delta = offset_vec(k); % positive offset
    
    % build the unit-max shape for THIS offset
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
    
    % scale to physical B1
    b1_single = B1_max * pulse_shape(:); % N×1, µT
    
    % distribute into the three-band table
    b1pulse = zeros(npoints, 3); % N×3
    half = b1_single / sqrt(2);
    b1pulse(:,1) = half; % –Δ
    b1pulse(:,3) = half; % +Δ
    
    % set the Delta_Hz vector (always ±delta, 0)
    Delta_Hz = [-delta 0 +delta];
    
    % propagate once through Bloch-McConnell
    Mz_vec_2band(k) = ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
    
end

%% Plot both cases on the same diagram
figure('DockControls','on');
hold on;

% Plot 1-band case
plot(offset_vec/1e3, Mz_vec_1band, 'b-', 'LineWidth', 1.6, 'DisplayName', '1-band');

% Plot 2-band case
plot(offset_vec/1e3, Mz_vec_2band, 'r-', 'LineWidth', 1.6, 'DisplayName', '2-band');

% Formatting
xlabel('\Delta (kHz)');
ylabel('$|M_{z,\mathrm{free}}|$ after RF pulse', 'Interpreter','latex');
title('Simulated Z-spectrum Comparison (9.4 T)');
legend('show', 'Location', 'best');
grid on;
hold off;

%% Calculate and plot ihMT difference
figure('DockControls','on');
ihMT_diff = Mz_vec_1band - Mz_vec_2band; % ihMT contrast
plot(offset_vec/1e3, ihMT_diff, 'g-', 'LineWidth', 1.6);
xlabel('\Delta (kHz)');
ylabel('ihMT contrast (1-band - 2-band)');
title('ihMT Contrast (9.4 T)');
grid on;

%% Display summary statistics
fprintf('\nSimulation completed!\n');
fprintf('Frequency range: %.1f Hz to %.1f kHz\n', offset_vec(1), offset_vec(end)/1e3);
fprintf('1-band Mz range: %.4f to %.4f\n', min(Mz_vec_1band), max(Mz_vec_1band));
fprintf('2-band Mz range: %.4f to %.4f\n', min(Mz_vec_2band), max(Mz_vec_2band));
fprintf('Maximum ihMT contrast: %.4f at %.1f kHz\n', max(ihMT_diff), offset_vec(ihMT_diff == max(ihMT_diff))/1e3);