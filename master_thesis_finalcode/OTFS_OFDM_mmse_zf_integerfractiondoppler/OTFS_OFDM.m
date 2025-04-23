%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%

clc
clear all
close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 4;
% number of subcarriers
M = 8;
% size of constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

SNR_dB = 0:2:20;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 10^4;
err_ber = zeros(length(SNR_dB),6);

for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        x = qammod(data_temp,M_mod,'gray');
        x = reshape(x,M,N);
        
        %% OTFS modulation%%%%
        s1 = OTFS_modulation(N,M,x);

        %% OTFS channel generation%%%%
        [taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen(N,M);
        [taps0,delay_taps0,Doppler_taps0,chan_coef0] = OTFS_channel_gen_fractionaldoppler(N,M);        
        
        %% OTFS channel output%%%%%
        r1 = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),s1);
        r2 = OTFS_channel_output(N,M,taps0,delay_taps0,Doppler_taps0,chan_coef0,sigma_2(iesn0),s1);        
        
        %% OTFS demodulation%%%%
        y1 = OTFS_demodulation(N,M,r1);
        y2 = OTFS_demodulation(N,M,r2);

        %% message passing detector%%%%
        
        x_est1 = OTFS_zf_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,y1);
        x_est2 = OTFS_mmse_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y1);
        x_est3 = OTFS_zf_detector(N,M,taps0,delay_taps0,Doppler_taps0,chan_coef0,y2);       
        x_est4 = OTFS_mmse_detector(N,M,taps0,delay_taps0,Doppler_taps0,chan_coef0,sigma_2(iesn0),y2); 
        
        %x_est5 = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y1); 
        %x_est6 = OFDM_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y3); 

        %% output bits and errors count%%%%%
        
        data_demapping1 = qamdemod(x_est1,M_mod,'gray');
        data_demapping2 = qamdemod(x_est2,M_mod,'gray');
        data_demapping3 = qamdemod(x_est3,M_mod,'gray');  
        data_demapping4 = qamdemod(x_est4,M_mod,'gray');
        
        %data_demapping5 = qamdemod(x_est5,M_mod,'gray');
        %data_demapping6 = qamdemod(x_est6,M_mod,'gray');

        
        data_info_est1 = reshape(de2bi(data_demapping1,M_bits),N_bits_perfram,1);
        data_info_est2 = reshape(de2bi(data_demapping2,M_bits),N_bits_perfram,1);
        data_info_est3 = reshape(de2bi(data_demapping3,M_bits),N_bits_perfram,1);
        data_info_est4 = reshape(de2bi(data_demapping4,M_bits),N_bits_perfram,1);
        
        %data_info_est5 = reshape(de2bi(data_demapping5,M_bits),N_bits_perfram,1); 
        %data_info_est6 = reshape(de2bi(data_demapping6,M_bits),N_bits_perfram,1);

        
        errors1 = sum(xor(data_info_est1,data_info_bit));
        errors2 = sum(xor(data_info_est2,data_info_bit));
        errors3 = sum(xor(data_info_est3,data_info_bit));        
        errors4 = sum(xor(data_info_est4,data_info_bit));
        
        %errors5 = sum(xor(data_info_est5,data_info_bit));   
        %errors6 = sum(xor(data_info_est6,data_info_bit)); 

        
        err_ber(iesn0,1) = errors1 + err_ber(iesn0,1);
        err_ber(iesn0,2) = errors2 + err_ber(iesn0,2);
        err_ber(iesn0,3) = errors3 + err_ber(iesn0,3); 
        err_ber(iesn0,4) = errors4 + err_ber(iesn0,4); 
        
        %err_ber(iesn0,5) = errors5 + err_ber(iesn0,5);  
        %err_ber(iesn0,6) = errors6 + err_ber(iesn0,6);  
        %ifram
    end
end


err_ber_fram1 = err_ber(:,1)/N_bits_perfram/N_fram;
err_ber_fram2 = err_ber(:,2)/N_bits_perfram/N_fram;
err_ber_fram3 = err_ber(:,3)/N_bits_perfram/N_fram;
err_ber_fram4 = err_ber(:,4)/N_bits_perfram/N_fram;

%err_ber_fram5 = err_ber(:,5)/N_bits_perfram/N_fram;
%err_ber_fram6 = err_ber(:,6)/N_bits_perfram/N_fram;



semilogy(SNR_dB, err_ber_fram1,'r-o','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram2,'g--o','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram3,'b-x','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram4,'k--x','LineWidth',2);
hold on

%semilogy(SNR_dB, err_ber_fram5,'k-s','LineWidth',2);
%hold on
%semilogy(SNR_dB, err_ber_fram6,'g-x','LineWidth',2);
%hold on

%legend('OTFS-ZF','OFDM-ZF','OTFS-MMSE','OFDM-MMSE','OTFS-MP')
%legend('OTFS-ZF','OTFS-MMSE','OTFS-MP')
%legend('MMSE-OTFS','MMSE-OFDM')
legend('OTFS-ZF-integerdoppler','OTFS-MMSE-integerdoppler','OTFS-ZF-fractionaldoppler','OTFS-MMSE-fractionaldoppler')
%title(sprintf('OTFS vs OFDM with MMSE / ZF detectors'));
title(sprintf('OTFS  with  ZF / MMSE detectors with integer / fractional doppler shift'));
ylabel('BER'); xlabel('SNR in dB');grid on

toc
