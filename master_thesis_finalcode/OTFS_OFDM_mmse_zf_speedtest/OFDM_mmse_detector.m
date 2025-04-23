%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function x_est = OFDM_mmse_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma,y)

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
        %FN(i,j) = exp(-1i*2*(pi)*(i-1)*(j-1)/N);
    end    
end

%Fn = reshape(FN,1,N*N);

IN = eye(N,N);

Heff = zeros(N*M,N*M);

for i= 1:taps
    dmatrix = zeros(N*M,N*M);
    for  j= 1:N*M
        dmatrix(j,j) = exp(1i*2*(pi)*Doppler_taps(i)*(j-1)/M/N);
    end
    %{
    %Heff = Heff + chan_coef(i)*(pmatrix^delay_taps(i))*dmatrix;
    %Heff = Heff + chan_coef(i)*FM*(pmatrix^delay_taps(i))*dmatrix*FM';

    %P = kron(FN,IN)*(pmatrix^delay_taps(i))*kron(FN',IN);
    %Q = kron(FN,IN)*dmatrix*kron(FN',IN);
    %Heff = Heff + chan_coef(i)*P*Q;
    %tmp = Heff + chan_coef(i)*kron(FM,IN)*(pmatrix^delay_taps(i))*dmatrix*kron(FM',IN);
    %tmp = Heff + chan_coef(i)*(pmatrix^delay_taps(i))*dmatrix;
    tmp = pmatrix^delay_taps(i);
    *Fn.'*(Fn').'
    tmp1 = diag(Fn);
    tmp2 = diag(Fn');
    tmp3 = diag(Fn)*diag(Fn');
    tmp4 = kron(IN',FN)*kron(IN,FN);
    P = kron(IN,FN)*(pmatrix^delay_taps(i))*kron(IN',FN);
    Q = kron(IN,FN)*dmatrix*kron(IN',FN); 
    Heff = Heff + chan_coef(i)*P*Q;
    %}

    P = kron(IN,FM)*(pmatrix^delay_taps(i))*kron(IN,FM');
    Q = kron(IN,FM)*dmatrix*kron(IN,FM'); 
    Heff = Heff + chan_coef(i)*P*Q;

end

I0 = eye(N*M,N*M);

%W=(Heff'*Heff+I0*sigma)\Heff';
W = pinv(Heff'*Heff+I0*sigma)*Heff';

x_est0 = W*yv;
x_est = reshape(x_est0,N,M);

end
