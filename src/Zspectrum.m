% ========================================================================
%  Z-spectrum simulation
%  –6 kHz … +6 kHz   •   63 sample points
% ========================================================================

%% pulse parameters
pulse_duration = 5; % 5 s saturation pulse
npoints        = 100000; % samples in the shape
B1_max         = 6.9; % µT (peak B1)
nband          = '1band';
shape          = 'square';
dt = pulse_duration/npoints; %dwell time

%% frequency offsets to sweep 
offset_vec = linspace(-6e3, 6e3, 63);      % Hz
Mz_vec     = zeros(size(offset_vec));          % pre-allocate result

%%  tissue parameters 
tissuepars           = init_tissue('hc');
tissuepars.lineshape = 'SL';


%%  sweep over offsets 
for k = 1:numel(offset_vec)

    delta    = abs(offset_vec(k));             % positive magnitude
    signOff  = sign(offset_vec(k));            % +1 or –1 (for 1-band)

    % build the unit-max shape for THIS offset 
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);

    % scale to physical B1 
    b1_single = B1_max * pulse_shape(:);       % N×1, µT

    % distribute into the three-band table 
    b1pulse = zeros(npoints, 3);               % N×3

    switch nband
        case '1band'                           % single side-band
            if signOff >= 0
                b1pulse(:,3) = b1_single;      % +Δ column only
            else
                b1pulse(:,1) = b1_single;      % –Δ column only
            end

        case '2band'                           % symmetrical ±Δ
            half = b1_single / sqrt(2);        
            b1pulse(:,1) = half;               % –Δ
            b1pulse(:,3) = half;               % +Δ
    end

    % set the Delta_Hz vector (always ±delta, 0)
    Delta_Hz = [-delta 0 +delta];

    % propagate once through Bloch-McConnell
    Mz_vec(k) = simplified_ssSPGR_ihMT_integrate(b1pulse, dt, Delta_Hz, tissuepars);
end

%% plot the Z-spectrum 
figure('DockControls','on');                   
plot(offset_vec/1e3, Mz_vec, 'LineWidth', 1.6);
xlabel('\Delta  (kHz)');
ylabel('$|M_{z,\mathrm{free}}|$ after RF pulse', 'Interpreter','latex');
title('Simulated Z-spectrum (9.4 T)');
grid on