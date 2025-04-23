%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function s = OFDM_modulation(N,M,x)
%% OTFS Modulation: 1. ISFFT, 2. Heisenberg transform
%X = fft(ifft(x).').'/sqrt(M/N); %%%ISFFT
s_mat = ifft(x)*sqrt(M); % Heisenberg transform
s = reshape(s_mat,N*M,1);
%s = reshape(x,N*M,1);
end
