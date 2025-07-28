%% Steady-State Test
% Test if given pulse duration reaches steady-state

clear; clc; close all;

%% Same parameters as z-spectra simulation
pulse_duration = 5; % saturation pulse duration
npoints = 100000; % samples in the shape
B1_max = 3.45*sqrt(2); % µT (peak B1)
nband = '2band';
dt = pulse_duration/npoints; %dwell time
shape = 'square';

%% Tissue parameters
tissuepars = init_tissue('hc');
tissuepars.lineshape = 'SL';

%% Test at a few different offsets
test_offsets = [6000]; % Hz

for i = 1:length(test_offsets)
    delta = test_offsets(i);
    fprintf('Testing at %.0f Hz offset:\n', delta);
    
    % Build pulse
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
    b1_single = B1_max * pulse_shape(:);
    b1pulse = zeros(npoints, 3);%create 3-band array
    switch nband
        case '1band'
            % For 1-band: put pulse in positive frequency band only
            signOff = 1; % assuming positive offset for this test
            if signOff >= 0
                b1pulse(:,3) = b1_single; % +Δ column only
            else
                b1pulse(:,1) = b1_single; % -Δ column only  
            end
            
        case '2band'
            half = b1_single / sqrt(2);        
            b1pulse(:,1) = half;               % –Δ
            b1pulse(:,3) = half;               % +Δ
            
    end
    
    Delta_Hz = [-delta 0 +delta];
    
    % Call your function with tracking and plotting
    [Mz_free, Mz_evolution] = steadystate_confirm_integrate(b1pulse, dt, Delta_Hz, tissuepars, 'track_convergence', 'plot_convergence');
    
    % Add title to distinguish plots
    title(sprintf('Steady-State Test at %.0f Hz', delta));
    
    fprintf('\n');
end

fprintf('All tests completed!\n');