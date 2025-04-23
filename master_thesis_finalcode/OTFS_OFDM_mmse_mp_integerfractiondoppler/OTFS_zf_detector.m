%Owner is Zuoyuan Fu, from Communication Systems, KTH.
%%
function x_est = OTFS_zf_detector(N,M,taps,delay_taps,Doppler_taps,chan_coef,y)

yv = reshape(y,N*M,1);
%yv = y(:);

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

%W=(Heff'*Heff+I0*sigma)\Heff';
W = pinv(Heff'*Heff)*Heff';

x_est = W*yv;
%x_est = reshape(x_est0,N,M);

end
