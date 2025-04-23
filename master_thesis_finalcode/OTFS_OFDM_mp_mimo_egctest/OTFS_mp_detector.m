

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
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%    - Latest version of this code may be downloaded from: https://ecse.monash.edu/staff/eviterbo/
%    - Freely distributed for educational and research purposes
%%
function x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,Heff,sigma_2,y,ni)

yv = reshape(y,N*M,1);
n_ite = 20;
delta_fra = 0.6;

alphabet = qammod(0:M_mod-1,M_mod,'gray');

mean_int = zeros(N*M,taps*(2*ni+1));

var_int = zeros(N*M,taps*(2*ni+1));

p_map = ones(N*M,taps*(2*ni+1),M_mod)*(1/M_mod);

conv_rate_prev = -0.1;

for ite=1:n_ite
    %% Update mean and var
    for ele1=1:1:N
        for ele2=1:1:M

            mean_int_hat = zeros(taps*(2*ni+1),1);
            var_int_hat = zeros(taps*(2*ni+1),1);

            mean_int_sum = 0;
            var_int_sum = 0;

            for tap_no=1:taps
                             
                %n = mod(round(ele1-Doppler_taps(tap_no)),N);

                for k = -ni:1:ni
                
                    m = mod(ele2-delay_taps(tap_no)-1,M)+1;
                    n = mod(round(ele1-Doppler_taps(tap_no))+k-1,N)+1;
    
                    %m = mod(round(ele2-Doppler_taps(tap_no)-1),M)+1;
                    %n = mod(ele1-delay_taps(tap_no)-1,N)+1;
    
                    new_chan = Heff(M*(ele1-1)+ele2 , M*(n-1)+m);

                    mean_tmp = 0;
                    var_tmp  = 0;
    
                    for i2=1:1:M_mod
                        mean_tmp = mean_tmp + p_map(M*(ele1-1)+ele2,(tap_no-1)*(2*ni+1)+k+ni+1,i2) * alphabet(i2);
                        var_tmp = var_tmp + p_map(M*(ele1-1)+ele2,(tap_no-1)*(2*ni+1)+k+ni+1,i2) * abs(alphabet(i2))^2;
                    end

                    mean_int_sum = mean_int_sum + mean_tmp * new_chan;
                    var_int_sum = var_int_sum + var_tmp* abs(new_chan)^2 - abs(mean_tmp* new_chan)^2;
                    % var_int_sum = var_int_sum + var_tmp* abs(new_chan)^2 - abs(mean_tmp)^2;

                    
                    mean_int_hat((tap_no-1)*(2*ni+1)+k+ni+1) = mean_tmp * new_chan;
                    var_int_hat((tap_no-1)*(2*ni+1)+k+ni+1) = var_tmp* abs(new_chan)^2 - abs(mean_tmp* new_chan)^2;
                       

                end               

            end
            
            %mean_int_sum = sum(mean_int_hat);
            var_int_sum = var_int_sum+(sigma_2);
            
            for tap_no=1:taps
                for k = -ni:1:ni
                    mean_int(M*(ele1-1)+ele2,(tap_no-1)*(2*ni+1)+k+ni+1) = mean_int_sum - mean_int_hat((tap_no-1)*(2*ni+1)+k+ni+1);
                    var_int(M*(ele1-1)+ele2,(tap_no-1)*(2*ni+1)+k+ni+1) = var_int_sum - var_int_hat((tap_no-1)*(2*ni+1)+k+ni+1);
                end
            end
            
        end
    end
    %% Update probabilities
    sum_prob_comp = zeros(N*M,M_mod);
    dum_eff_ele1 = zeros(taps*(2*ni+1),1);
    dum_eff_ele2 = zeros(taps*(2*ni+1),1);

    for ele1=1:1:N
        for ele2=1:1:M

            dum_sum_prob = zeros(M_mod,1);
            log_te_var = zeros(taps*(2*ni+1),M_mod);

            for tap_no=1:taps

                for k = -ni:1:ni

                %eff_ele1 = mod(ele1+delay_taps(tap_no)-1,N)+1;
                %eff_ele2 = mod(round(ele2+Doppler_taps(tap_no)-1),M)+1;  

                    eff_ele1 = mod(round(ele1+Doppler_taps(tap_no)-k-1),N)+1;
                    eff_ele2 = mod(ele2+delay_taps(tap_no)-1,M)+1;
    
                    new_chan = Heff(M*(eff_ele1-1)+eff_ele2,M*(ele1-1)+ele2); 
    
                    %new_chan = Heff(M*(ele1-1)+ele2 , M*n + m + 1);
                    
                    dum_eff_ele1((tap_no-1)*(2*ni+1)+k+ni+1) = eff_ele1;
                    dum_eff_ele2((tap_no-1)*(2*ni+1)+k+ni+1) = eff_ele2;
                    for i2=1:1:M_mod
                        dum_sum_prob(i2) = abs(yv(M*(eff_ele1-1)+eff_ele2) - mean_int(M*(eff_ele1-1)+eff_ele2,(tap_no-1)*(2*ni+1)+k+ni+1) - new_chan * alphabet(i2))^2;
                        dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(M*(eff_ele1-1)+eff_ele2,(tap_no-1)*(2*ni+1)+k+ni+1));
                    end
                    dum_sum = dum_sum_prob - max(dum_sum_prob);
                    dum1 = sum(exp(dum_sum));
                    log_te_var((tap_no-1)*(2*ni+1)+k+ni+1,:) = dum_sum - log(dum1);

                end
                
            end

            for i2=1:1:M_mod
                ln_qi(i2) = sum(log_te_var(:,i2));
            end
            dum_sum = exp(ln_qi - max(ln_qi));
            dum1 = sum(dum_sum);
            sum_prob_comp(M*(ele1-1)+ele2,:) = dum_sum/dum1;

            for tap_no=1:1:taps
                for k = -ni:1:ni

                    eff_ele1 = dum_eff_ele1((tap_no-1)*(2*ni+1)+k+ni+1);
                    eff_ele2 = dum_eff_ele2((tap_no-1)*(2*ni+1)+k+ni+1);
                    
                    dum_sum = log_te_var((tap_no-1)*(2*ni+1)+k+ni+1,:);
                    ln_qi_loc = ln_qi - dum_sum;
                    dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                    dum1 = sum(dum_sum);
                    p_map(M*(eff_ele1-1)+eff_ele2,(tap_no-1)*(2*ni+1)+k+ni+1,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(M*(eff_ele1-1)+eff_ele2,(tap_no-1)*(2*ni+1)+k+ni+1,:),1,M_mod);

                end

            end
            
        end
    end
    conv_rate =  sum(max(sum_prob_comp,[],2)>0.99)/(N*M);
    if conv_rate==1
        sum_prob_fin = sum_prob_comp;
        break;
    elseif conv_rate > conv_rate_prev
        conv_rate_prev = conv_rate;
        sum_prob_fin = sum_prob_comp;
    elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
        break;
    end
end
x_est = zeros(M,N);
for ele1=1:1:N
    for ele2=1:1:M
        [~,pos] = max(sum_prob_fin(M*(ele1-1)+ele2,:));
        x_est(ele2,ele1) = alphabet(pos);
    end
end
end