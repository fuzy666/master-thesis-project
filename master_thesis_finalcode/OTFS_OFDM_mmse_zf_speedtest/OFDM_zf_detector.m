%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function x_est = OFDM_zf_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,y)

yv = reshape(y,N*M,1);
%yv = y(:);

pmatrix = zeros(N*M,N*M);
pmatrix(1,N*M) = 1;
for i = 2:N*M
    pmatrix(i,i-1) = 1;
end


FM = zeros(M,M);
for i = 1:M
    for j = 1:M
        FM(i,j) = exp(-1i*2*(pi)*(i-1)*(j-1)/M)/sqrt(M);
    end    
end

IN = eye(N,N);

Heff = zeros(N*M,N*M);

for i= 1:taps
    dmatrix = zeros(N*M,N*M);
    for  j= 1:N*M
        dmatrix(j,j) = exp(1i*2*(pi)*Doppler_taps(i)*(j-1)/M/N);
    end

    P = kron(IN,FM)*(pmatrix^delay_taps(i))*kron(IN,FM');
    Q = kron(IN,FM)*dmatrix*kron(IN,FM'); 
    Heff = Heff + chan_coef(i)*P*Q;

end

%W=(Heff'*Heff+I0*sigma)\Heff';
W = pinv(Heff'*Heff)*Heff';

x_est0 = W*yv;
x_est = reshape(x_est0,N,M);

end
