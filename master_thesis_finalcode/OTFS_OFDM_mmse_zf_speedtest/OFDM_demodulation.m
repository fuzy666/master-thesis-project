%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function y = OFDM_demodulation(N,M,r)
%% OTFS demodulation: 1. Wiegner transform, 2. SFFT
r_mat = reshape(r,M,N);
y = fft(r_mat)/sqrt(M); % Wigner transform
%y = r_mat;
end
