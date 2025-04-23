%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,y)

yv = reshape(y,N*M,1);
n_ite = 50;
delta_fra = 0.6;
ni = 1;

alphabet = qammod(0:M_mod-1,M_mod,'gray');

mean_int = zeros(N*M,taps);
var_int = zeros(N*M,taps);
p_map = ones(N*M,taps,M_mod)*(1/M_mod);

conv_rate_prev = -0.1;

pmatrix = zeros(N*M,N*M);
pmatrix(1,N*M) = 1;
for i = 2:N*M
    pmatrix(i,i-1) = 1;
end

FN = zeros(N,N);
for i = 1:N
    for j = 1:N
        FN(i,j) = exp(-1i*2*(pi)*(i-1)*(j-1)/N)/sqrt(N);
    end    
end

IM = eye(M,M);

Heff = zeros(N*M,N*M);

for i= 1:taps

    dmatrix = zeros(N*M,N*M);
    for j = 1:N*M
        dmatrix(j,j) = exp(1i*2*(pi)*Doppler_taps(i)*(j-1)/M/N);
    end

    P = kron(FN,IM)*(pmatrix^delay_taps(i))*kron(FN',IM);
    Q = kron(FN,IM)*dmatrix*kron(FN',IM);
    Heff = Heff + chan_coef(i)*P*Q;

end

for ite=1:n_ite
    %% Update mean and var
    for ele1=1:1:N
        for ele2=1:1:M

            mean_int_hat = zeros(taps,1);
            var_int_hat = zeros(taps,1);

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
                        mean_tmp = mean_tmp + p_map(M*(ele1-1)+ele2,tap_no,i2) * alphabet(i2);
                        var_tmp = var_tmp + p_map(M*(ele1-1)+ele2,tap_no,i2) * abs(alphabet(i2))^2;
                    end

                    mean_int_sum = mean_int_sum + mean_tmp * new_chan;
                    var_int_sum = var_int_sum + var_tmp* abs(new_chan)^2 - abs(mean_tmp* new_chan)^2;
                    % var_int_sum = var_int_sum + var_tmp* abs(new_chan)^2 - abs(mean_tmp)^2;

                    if k == 0
                        mean_int_hat(tap_no) = mean_tmp * new_chan;
                        var_int_hat(tap_no) = var_tmp* abs(new_chan)^2 - abs(mean_tmp* new_chan)^2;
                    end   

                end               

            end
            
            %mean_int_sum = sum(mean_int_hat);
            var_int_sum = var_int_sum+(sigma_2);
            
            for tap_no=1:taps
                mean_int(M*(ele1-1)+ele2,tap_no) = mean_int_sum - mean_int_hat(tap_no);
                var_int(M*(ele1-1)+ele2,tap_no) = var_int_sum - var_int_hat(tap_no);
            end
            
        end
    end
    %% Update probabilities
    sum_prob_comp = zeros(N*M,M_mod);
    dum_eff_ele1 = zeros(taps,1);
    dum_eff_ele2 = zeros(taps,1);

    for ele1=1:1:N
        for ele2=1:1:M

            dum_sum_prob = zeros(M_mod,1);
            log_te_var = zeros(taps,M_mod);

            for tap_no=1:taps

                %eff_ele1 = mod(ele1+delay_taps(tap_no)-1,N)+1;
                %eff_ele2 = mod(round(ele2+Doppler_taps(tap_no)-1),M)+1;  

                eff_ele1 = mod(round(ele1+Doppler_taps(tap_no)-1),N)+1;
                eff_ele2 = mod(ele2+delay_taps(tap_no)-1,M)+1;

                new_chan = Heff(M*(eff_ele1-1)+eff_ele2,M*(ele1-1)+ele2); 

                %new_chan = Heff(M*(ele1-1)+ele2 , M*n + m + 1);
                
                dum_eff_ele1(tap_no) = eff_ele1;
                dum_eff_ele2(tap_no) = eff_ele2;
                for i2=1:1:M_mod
                    dum_sum_prob(i2) = abs(yv(M*(eff_ele1-1)+eff_ele2) - mean_int(M*(eff_ele1-1)+eff_ele2,tap_no) - new_chan * alphabet(i2))^2;
                    dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(M*(eff_ele1-1)+eff_ele2,tap_no));
                end
                dum_sum = dum_sum_prob - max(dum_sum_prob);
                dum1 = sum(exp(dum_sum));
                log_te_var(tap_no,:) = dum_sum - log(dum1);
                
            end

            for i2=1:1:M_mod
                ln_qi(i2) = sum(log_te_var(:,i2));
            end
            dum_sum = exp(ln_qi - max(ln_qi));
            dum1 = sum(dum_sum);
            sum_prob_comp(M*(ele1-1)+ele2,:) = dum_sum/dum1;

            for tap_no=1:1:taps
                eff_ele1 = dum_eff_ele1(tap_no);
                eff_ele2 = dum_eff_ele2(tap_no);
                
                dum_sum = log_te_var(tap_no,:);
                ln_qi_loc = ln_qi - dum_sum;
                dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                dum1 = sum(dum_sum);
                p_map(M*(eff_ele1-1)+eff_ele2,tap_no,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(M*(eff_ele1-1)+eff_ele2,tap_no,:),1,M_mod);
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
