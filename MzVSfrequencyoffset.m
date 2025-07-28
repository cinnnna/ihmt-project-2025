%% ------------  Scanner / sequence parameters  --------------------------
pulse_duration = 3; % pulse duration (s)
b1max = 6.9; % µT
offset = 5500; % band offset (Hz)
nband = '2band';
npoints = 100000;

dt = pulse_duration/npoints; %dwell time

[pulse_shape] = te_gen_MB_pulse(pulse_duration,npoints,offset,nband,'square');


b1_single = b1max * pulse_shape(:);

% Make b1pulse into 3 columns form as the integral expects (keep column 2
% with zeros)
b1pulse = zeros(npoints, 3);  % N×3 as the integrator expects
switch nband
    case '1band'                             % single side-band
        if offset > 0
            b1pulse(:,3) = b1_single;        % +Δ only
        else
            b1pulse(:,1) = b1_single;        % –Δ only
        end

    case '2band'                             % symmetrical ±Δ
        half = b1_single / sqrt(2);          % split amplitude equally %divide by sqrt2
        b1pulse(:,1) = half;                 % –Δ column
        b1pulse(:,3) = half;                 % +Δ column

    otherwise
        error('nband must be ''1band'' or ''2band''.');
end


%% ------------  Frequency-offset sweep  ----------------------------------
Delta_vec = linspace(-55e2, 55e2, 63);   % Hz
Mz_vec    = zeros(size(Delta_vec));

%% ------------  Tissue parameters  --------------------------------------
tissuepars = init_tissue('BSA');   
fprintf('tissuepars.semi.T2 = %e\n', tissuepars.semi.T2);  
fprintf('All semi fields:\n');
disp(tissuepars.semi);
tissuepars.lineshape = 'Gaussian';      % 'SL' or 'Gaussian'


%% ------------  Loop over offsets  --------------------------------------
for k = 1:numel(Delta_vec)
    Delta_Hz = [-Delta_vec(k), 0, +Delta_vec(k)];   % three-band vector
    Mz_vec(k) = ssSPGR_ihMT_integrate( ...
                     b1pulse, dt, Delta_Hz, tissuepars);
end

%% ------------  Plot  ----------------------------------------------------
figure;
plot(Delta_vec/1e3, Mz_vec, 'LineWidth', 1.6);
xlabel('\Delta  (kHz)');
ylabel({'Mz_free', 'after RF pulse'}, ...
       'Interpreter','tex');   
title('Free-water longitudinal magnetisation vs. frequency offset');
grid on;
