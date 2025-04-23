%
% Copyright (c) 2018, Raviteja Patchava, Yi Hong, and Emanuele Viterbo, Monash University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR

%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSE
% QUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%    - Latest version of this code may be downloaded from: https://ecse.monash.edu/staff/eviterbo/
%    - Freely distributed for educational and research purposes
%%

clc
clear all
close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 16;
% number of subcarriers
M = 4;
% size of constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

ni_0 = 0:1:floor((N-1)/2);

SNR_dB = 20:2:20;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 10^4;

%err_ber = zeros(length(SNR_dB),6);
err_ber = zeros(length(ni_0),1);

ni = 0;

%for iesn0 = 1:length(SNR_dB)
for iesn0 = 1:length(ni_0)

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
        
        %% OTFS channel output%%%%%
        r1 = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,s1);
        
        %% OTFS demodulation%%%%
        y1 = OTFS_demodulation(N,M,r1);
        
        %% message passing detector%%%%
        x_est1 = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,y1,ni); 

        %% output bits and errors count%%%%%
        data_demapping1 = qamdemod(x_est1,M_mod,'gray');

        data_info_est1 = reshape(de2bi(data_demapping1,M_bits),N_bits_perfram,1); 

        errors1 = sum(xor(data_info_est1,data_info_bit));   

        err_ber(iesn0,1) = errors1 + err_ber(iesn0,1);  

    end

    ni = ni + 1;
end

err_ber_fram1 = err_ber(:,1)/N_bits_perfram/N_fram;

%semilogy(SNR_dB, err_ber_fram5,'k-s','LineWidth',2);
semilogy(ni_0, err_ber_fram1,'k-s','LineWidth',2);
hold on
%semilogy(SNR_dB, err_ber_fram6,'g-x','LineWidth',2);
%hold on

legend('OTFS-MP-ni-test')
title(sprintf('OTFS  with MP detectors with diff doppler interf'));
ylabel('BER'); xlabel('ni');grid on

toc