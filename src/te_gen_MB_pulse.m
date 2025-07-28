function [pulse_shape] = te_gen_MB_pulse(pulse_duration,npoints,offset,nband,shape)
%%% function based on publication
%%% [pulse] = te_gen_MB_pulse(pulse_duration,npoints,offset,nband)
%%%   	arguments:
% % 	-  pulse_duration:   pulse duration (s)
% % 	-  offset:	         offset frequency for off-resonant bands (Hz) 
% % 	-  nband:		     string argument for number of bands. Options are '1band', or '2band'.
%%% Feb 2019
%%% April 2009 - output the TBP as well
%%% 
%%
gausswin_alpha = 3; % shape of basic pulse, width of gaussian
dt = pulse_duration/npoints;
t=linspace(dt,pulse_duration,npoints);

%% basic shape

switch shape
    case 'gaussian'
        pulse = gausswin(npoints,gausswin_alpha);
    case 'square'
        pulse = ones(npoints,1);
end

%% Make pulse 

%%% Modulation functions
%switch nband
    %case '1band'
        %wt = exp(1i*2*pi*offset*t);
        
    %case '2band'
        %wt = cos(2*pi*offset*t);   
%end

%%% Make modulated pulse
%pulse_shape = pulse(:) .* wt(:);
pulse_shape = pulse (:);

%plot(t,pulse_shape);

end