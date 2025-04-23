%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%

clc
clear all
close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 8;
% number of subcarriers
M = 16;
% size of constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

SNR_dB = 20:2:20;
delta = 0.1:0.1:0.9;

SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 10^4;
%err_ber = zeros(length(SNR_dB),1);
err_ber = zeros(length(delta),1);

%for iesn0 = 1:length(SNR_dB)
for k = 1:length(delta)
    for ifram = 1:N_fram
        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        x = qammod(data_temp,M_mod,'gray');
        x = reshape(x,N,M);
        
        %% OTFS modulation%%%%
        s = OTFS_modulation(N,M,x);
        
        %% OTFS channel generation%%%%
        [taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen(N,M);
        
        %% OTFS channel output%%%%%
        r = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,s);
        
        %% OTFS demodulation%%%%
        y = OTFS_demodulation(N,M,r);
        
        %% message passing detector%%%%
        x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,y,delta(k));
        
        %% output bits and errors count%%%%%
        data_demapping = qamdemod(x_est,M_mod,'gray');
        data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
        errors = sum(xor(data_info_est,data_info_bit));
        err_ber(k) = errors + err_ber(k);
        %err_ber(SNR_dB) = errors + err_ber(SNR_dB);
        %ifram
    end
end
err_ber_fram = err_ber/N_bits_perfram./N_fram;
%err_ber_fram(2) = 0.000246406250000000;
semilogy(delta, err_ber_fram,'-*','LineWidth',2);
title(sprintf('OTFS with MP detector ( SNR = 20 dB and iterations number = 80 )'))
ylabel('BER'); xlabel('\Delta');grid on

toc
