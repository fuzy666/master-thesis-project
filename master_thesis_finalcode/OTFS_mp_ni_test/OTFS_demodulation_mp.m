%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function y = OTFS_demodulation_mp(N,M,r)
%% OTFS demodulation: 1. Wiegner transform, 2. SFFT

%r_mat = reshape(r,M,N);

r_mat = reshape(r,N,M);
Y = fft(r_mat)/sqrt(M); % Wigner transform
Y = Y.';
%y = ifft(fft(Y).')/sqrt(N/M); % SFFT

y = ifft(fft(Y).').'/sqrt(N/M); % SFFT

end
