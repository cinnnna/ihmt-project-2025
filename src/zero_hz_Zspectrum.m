% ========================================================================
% Z-spectrum simulation - MODIFIED FOR 2-BAND SYSTEM
% –6 kHz … +6 kHz • 63 sample points
% ========================================================================
%% pulse parameters
pulse_duration = 5; % 5 s saturation pulse
npoints = 100000; % samples in the shape
B1_max = 6.9*sqrt(2); % µT (peak B1)
nband = '2band';
shape = 'square';
dt = pulse_duration/npoints; %dwell time

%% frequency offsets to sweep
offset_vec = linspace(0, 6e3, 31); % Hz
Mz_vec = zeros(size(offset_vec)); % pre-allocate result

%% tissue parameters
tissuepars = init_tissue('BSA');
tissuepars.lineshape = 'Gaussian';

%% sweep over offsets
for k = 1:numel(offset_vec)
    delta = abs(offset_vec(k)); % positive magnitude

    % build the unit-max shape for THIS offset
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
    
    b1_band = B1_max * pulse_shape(:); % N×1, µT

    % fprintf('Pulse shape (first value): %.3f + %.3fi\n', real(b1_band(1)), imag(b1_band(1)));
    % fprintf('Pulse shape (second value): %.3f + %.3fi\n', real(b1_band(2)), imag(b1_band(2)));

    % propagate once through Bloch-McConnell (MODIFIED FUNCTION NAME)
    Mz_vec(k) = new_Dualcase_ssSPGR_ihMT_integrate(b1_band, dt, delta, tissuepars, nband);
end


%% plot the Z-spectrum
figure('DockControls','on');
plot(offset_vec/1e3, Mz_vec, 'LineWidth', 1.6);
xlabel('\Delta (kHz)');
ylabel('$|M_{z,\mathrm{free}}|$ after RF pulse', 'Interpreter','latex');
title('Simulated Z-spectrum (9.4 T)');
grid on;