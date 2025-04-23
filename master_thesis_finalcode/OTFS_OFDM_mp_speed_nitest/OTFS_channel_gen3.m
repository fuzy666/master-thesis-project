%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function [taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen3()
%% Channel for testing%%%%%
%channel with 4 taps of uniform power%%% 
taps = 4;
delay_taps = [0 1 2 3];
Doppler_taps = [137.7 139.3 141.7 142.9];
pow_prof = (1/taps) * (ones(1,taps));
chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
%%%%%%%%%%%%%%%%%%%%

end
