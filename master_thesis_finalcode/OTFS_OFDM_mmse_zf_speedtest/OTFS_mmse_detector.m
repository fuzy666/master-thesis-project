function x_est = OTFS_mmse_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma,y)

yv = reshape(y,N*M,1);
%yv = y(:);

pmatrix = zeros(N*M,N*M);
pmatrix(1,N*M) = 1;
for i = 2:N*M
    pmatrix(i,i-1) = 1;
end

%dmatrix = zeros(N*M,N*M);
%for i = 1:N*M
%    dmatrix(i,i) = exp(1i*2*(pi/M)*(i-1)/N);
%end

FN = zeros(N,N);
for i = 1:N
    for j = 1:N
        FN(i,j) = exp(-1i*2*(pi)*(i-1)*(j-1)/N)/sqrt(N);
        %FN(i,j) = exp(-1i*2*(pi)*(i-1)*(j-1)/N);
    end    
end

IM = eye(M,M);

%tmp = FN*FN';
Heff = zeros(N*M,N*M);
%tmp = zeros(N*M,N*M);
%tmp1 = zeros(N*M,N*M);
%tmp2 = zeros(N*M,N*M);

for i= 1:taps

    dmatrix = zeros(N*M,N*M);
    for j = 1:N*M
        dmatrix(j,j) = exp(1i*2*(pi)*Doppler_taps(i)*(j-1)/M/N);
    end

    P = kron(FN,IM)*(pmatrix^delay_taps(i))*kron(FN',IM);
    Q = kron(FN,IM)*dmatrix*kron(FN',IM);
    %T = P*Q/(N^2);
    Heff = Heff + chan_coef(i)*P*Q;

    %{
    T1 = kron(FN,IM);
    T2 = kron(FN',IM);
    T3 = kron(IM,FN);
    t1 = FN;
    t2 = FN';
    t3 = FN*FN';
    tmp = tmp + chan_coef(i)*kron(FN,IM)*(pmatrix^delay_taps(i))*(dmatrix^Doppler_taps(i))*kron(FN',IM);
    tmp = dmatrix^Doppler_taps(i);
    tmp1 = kron(FN,IM)*(pmatrix^delay_taps(i))*kron(FN',IM);
    tmp2 = kron(FN,IM)*(dmatrix^Doppler_taps(i))*kron(FN',IM);
    %}

end

I0 = eye(N*M,N*M);

%Heff0 = kron(FN,IM)*hmatrix*kron(FN',IM);

%W=(Heff'*Heff+I0*sigma)\Heff';
W = pinv(Heff'*Heff+I0*sigma)*Heff';

x_est0 = W*yv;
x_est = reshape(x_est0,N,M);

end