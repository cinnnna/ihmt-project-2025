% ========================================================================
% Z-spectrum simulation - MODIFIED FOR 2-BAND SYSTEM
% –6 kHz … +6 kHz • 63 sample points
% ========================================================================
%% pulse parameters
pulse_duration = 5; % 5 s saturation pulse
npoints = 100000; % samples in the shape
B1_max = 6.9; % µT (peak B1)
nband = '1band';
shape = 'square';
dt = pulse_duration/npoints; %dwell time

%% frequency offsets to sweep
offset_vec = linspace(-6e3, 6e3, 63); % Hz
Mz_vec = zeros(size(offset_vec)); % pre-allocate result

%% tissue parameters
tissuepars = init_tissue('hc');
tissuepars.lineshape = 'SL';

%% sweep over offsets
for k = 1:numel(offset_vec)
    delta = abs(offset_vec(k)); % positive magnitude
    signOff = sign(offset_vec(k)); % +1 or –1 (for 1-band)
    
    % build the unit-max shape for THIS offset
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
    
    % scale to physical B1
    b1_single = B1_max * pulse_shape(:); % N×1, µT
    
    fprintf('Pulse shape values:\n');
    disp(pulse_shape(:));
    % % distribute into the TWO-band table (MODIFIED)
    % b1pulse = zeros(npoints, 2); % N×2 (changed from N×3)
    % 
    % % Alternative approach - always put pulse at the actual offset
    % if delta == 0
    %     % On-resonance case
    %     b1pulse(:,1) = b1_single; % 0 Hz band
    %     b1pulse(:,2) = 0;
    %     Delta_Hz = [0, 0];
    % else
    %     % For any non-zero offset, put pulse at that exact frequency
    %     if signOff >= 0
    %         % Positive offset: [0, +delta] with pulse at +delta
    %         b1pulse(:,1) = 0;
    %         b1pulse(:,2) = b1_single;
    %         Delta_Hz = [0, delta];
    %     else
    %         % Negative offset: [-delta, 0] with pulse at -delta  
    %         b1pulse(:,1) = b1_single;
    %         b1pulse(:,2) = 0;
    %         Delta_Hz = [-delta, 0];
    %     end
    % end
    
    % propagate once through Bloch-McConnell (MODIFIED FUNCTION NAME)
    Mz_vec(k) = newDualcase_ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);

end

%% plot the Z-spectrum
figure('DockControls','on');
plot(offset_vec/1e3, Mz_vec, 'LineWidth', 1.6);
xlabel('\Delta (kHz)');
ylabel('$|M_{z,\mathrm{free}}|$ after RF pulse', 'Interpreter','latex');
title('Simulated Z-spectrum (9.4 T)');
grid on;