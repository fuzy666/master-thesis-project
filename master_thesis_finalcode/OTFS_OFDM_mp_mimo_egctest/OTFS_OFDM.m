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

ni=2;

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
        %s2 = OFDM_modulation(N,M,x);

        %% OTFS channel generation%%%%
        [taps1,delay_taps1,Doppler_taps1,chan_coef1] = OTFS_channel_gen1();
        [taps2,delay_taps2,Doppler_taps2,chan_coef2] = OTFS_channel_gen2();
        [taps3,delay_taps3,Doppler_taps3,chan_coef3] = OTFS_channel_gen3();
        
        %% OTFS channel output%%%%%
        r1 = OTFS_channel_output(N,M,taps1,delay_taps1,Doppler_taps1,chan_coef1,sigma_2(iesn0),s1);
        r2 = OTFS_channel_output(N,M,taps2,delay_taps2,Doppler_taps2,chan_coef2,sigma_2(iesn0),s1);
        r3 = OTFS_channel_output(N,M,taps3,delay_taps3,Doppler_taps3,chan_coef3,sigma_2(iesn0),s1);
        %r2 = OFDM_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),s2);
        
        %% OTFS demodulation%%%%
        y1 = OTFS_demodulation(N,M,r1);
        y2 = OTFS_demodulation(N,M,r2);
        y3 = OTFS_demodulation(N,M,r3);
        %y2 = OFDM_demodulation(N,M,r2);  
        
        %% message passing detector%%%%
        %{
        x_est1 = OTFS_zf_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,y1);
        x_est2 = OFDM_zf_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,y2);
        x_est3 = OTFS_mmse_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y1);       
        x_est4 = OFDM_mmse_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y2); 
        %}
        Heff1 = OTFS_H(N,M,taps1,delay_taps1,Doppler_taps1,chan_coef1);   %[taps,delay_taps,Doppler_taps,chan_coef]
        Heff2 = OTFS_H(N,M,taps2,delay_taps2,Doppler_taps2,chan_coef2);  %[taps,delay_taps,Doppler_taps,chan_coef]
        Heff3 = OTFS_H(N,M,taps3,delay_taps3,Doppler_taps3,chan_coef3);  %[taps,delay_taps,Doppler_taps,chan_coef]

        x_est1 = OTFS_mp_detector(N,M,M_mod,taps1,delay_taps1,Doppler_taps1,Heff1,sigma_2(iesn0),y1,ni); 
        x_est2 = OTFS_mp_detector(N,M,M_mod,taps2,delay_taps2,Doppler_taps2,Heff2,sigma_2(iesn0),y2,ni); 
        x_est3 = OTFS_mp_detector(N,M,M_mod,taps3,delay_taps3,Doppler_taps3,Heff3,sigma_2(iesn0),y3,ni); 
        x_est4 = OTFS_mp_detector1(N,M,M_mod,taps1+taps2+taps3,union(union(delay_taps1,delay_taps2),delay_taps3),union(union(Doppler_taps1,Doppler_taps2),Doppler_taps3),Heff1+Heff2+Heff3,sigma_2(iesn0),y1+y2+y3,ni);         
        %x_est6 = OFDM_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y3); 

        %% output bits and errors count%%%%%
        %{
        data_demapping1 = qamdemod(x_est1,M_mod,'gray');
        data_demapping2 = qamdemod(x_est2,M_mod,'gray');
        data_demapping3 = qamdemod(x_est3,M_mod,'gray');  
        data_demapping4 = qamdemod(x_est4,M_mod,'gray');
        %}
        data_demapping1 = qamdemod(x_est1,M_mod,'gray');
        data_demapping2 = qamdemod(x_est2,M_mod,'gray');
        data_demapping3 = qamdemod(x_est3,M_mod,'gray');
        data_demapping4 = qamdemod(x_est4,M_mod,'gray');
        %data_demapping6 = qamdemod(x_est6,M_mod,'gray');

        %{
        data_info_est1 = reshape(de2bi(data_demapping1,M_bits),N_bits_perfram,1);
        data_info_est2 = reshape(de2bi(data_demapping2,M_bits),N_bits_perfram,1);
        data_info_est3 = reshape(de2bi(data_demapping3,M_bits),N_bits_perfram,1);
        data_info_est4 = reshape(de2bi(data_demapping4,M_bits),N_bits_perfram,1);
        %}
        data_info_est1 = reshape(de2bi(data_demapping1,M_bits),N_bits_perfram,1); 
        data_info_est2 = reshape(de2bi(data_demapping2,M_bits),N_bits_perfram,1); 
        data_info_est3 = reshape(de2bi(data_demapping3,M_bits),N_bits_perfram,1); 
        data_info_est4 = reshape(de2bi(data_demapping4,M_bits),N_bits_perfram,1);         
        %data_info_est6 = reshape(de2bi(data_demapping6,M_bits),N_bits_perfram,1);

        %{
        errors1 = sum(xor(data_info_est1,data_info_bit));
        errors2 = sum(xor(data_info_est2,data_info_bit));
        errors3 = sum(xor(data_info_est3,data_info_bit));        
        errors4 = sum(xor(data_info_est4,data_info_bit));
        %}
        errors1 = sum(xor(data_info_est1,data_info_bit)); 
        errors2 = sum(xor(data_info_est2,data_info_bit));
        errors3 = sum(xor(data_info_est3,data_info_bit));  
        errors4 = sum(xor(data_info_est4,data_info_bit));
        %errors6 = sum(xor(data_info_est6,data_info_bit)); 

        %{
        err_ber(iesn0,1) = errors1 + err_ber(iesn0,1);
        err_ber(iesn0,2) = errors2 + err_ber(iesn0,2);
        err_ber(iesn0,3) = errors3 + err_ber(iesn0,3); 
        err_ber(iesn0,4) = errors4 + err_ber(iesn0,4); 
        %}
        err_ber(iesn0,1) = errors1 + err_ber(iesn0,1);  
        err_ber(iesn0,2) = errors2 + err_ber(iesn0,2); 
        err_ber(iesn0,3) = errors3 + err_ber(iesn0,3); 
        err_ber(iesn0,4) = errors4 + err_ber(iesn0,4); 
        %err_ber(iesn0,6) = errors6 + err_ber(iesn0,6);  
        %ifram
    end
end

%{
err_ber_fram1 = err_ber(:,1)/N_bits_perfram/N_fram;
err_ber_fram2 = err_ber(:,2)/N_bits_perfram/N_fram;
err_ber_fram3 = err_ber(:,3)/N_bits_perfram/N_fram;
err_ber_fram4 = err_ber(:,4)/N_bits_perfram/N_fram;
%}
err_ber_fram1 = err_ber(:,1)/N_bits_perfram/N_fram;
err_ber_fram2 = err_ber(:,2)/N_bits_perfram/N_fram;
err_ber_fram3 = err_ber(:,3)/N_bits_perfram/N_fram;
err_ber_fram4 = err_ber(:,4)/N_bits_perfram/N_fram;
%err_ber_fram6 = err_ber(:,6)/N_bits_perfram/N_fram;


%{
semilogy(SNR_dB, err_ber_fram1,'r-o','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram2,'g--o','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram3,'b-x','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram4,'k--x','LineWidth',2);
hold on
%}
semilogy(SNR_dB, err_ber_fram1,'k-x','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram2,'b-x','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram3,'r-x','LineWidth',2);
hold on
semilogy(SNR_dB, err_ber_fram4,'g-x','LineWidth',2);
hold on
%semilogy(SNR_dB, err_ber_fram6,'g-x','LineWidth',2);
%hold on

%legend('OTFS-ZF','OFDM-ZF','OTFS-MMSE','OFDM-MMSE','OTFS-MP')
%legend('OTFS-ZF','OTFS-MMSE','OTFS-MP')
%legend('MMSE-OTFS','MMSE-OFDM')
legend('OTFS-chan1','OTFS-chan2','OTFS-chan3','OTFS-multi-chans')
%title(sprintf('OTFS vs OFDM with MMSE / ZF detectors'));
title(sprintf('OTFS  with  single-chan/multi-chans'));
ylabel('BER'); xlabel('SNR in dB');grid on

toc
