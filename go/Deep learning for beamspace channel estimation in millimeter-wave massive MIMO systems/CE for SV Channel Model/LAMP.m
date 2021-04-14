function[xhat]=LAMP(y,A,type,snr_range)
% y:measurement
% A:sensing matrix

[M,N]=size(A);
T=8; % the number of layer
eta_l1=@(r,lam) (r./abs(r)).*max(bsxfun(@minus,abs(r),lam),0); %soft-thresholding shrinkage function

% xhat=zeros(N,L); %Initialization of signal estimation
% nmse=zeros(1,T);

if type==1 && snr_range==1
    load(['Trained_LAMP_for_SV_ULA_',num2str(N),num2str(M),'0to10dB.mat'])
end
if type==1 && snr_range==2
    load(['Trained_LAMP_for_SV_ULA_',num2str(N),num2str(M),'10to20dB.mat'])
end
if type==2 && snr_range==1
    load(['Trained_LAMP_for_SV_UPA_',num2str(N),num2str(M),'0to10dB.mat'])
end
if type==2 && snr_range==2
    load(['Trained_LAMP_for_SV_UPA_',num2str(N),num2str(M),'10to20dB.mat'])
end

v=y;  %Initialization of residual
theta=theta0(1);
B=B_0;
rvar = (1/M)*sum(abs(y).^2,1);
lam=theta*sqrt(rvar);
lam=max(abs(lam),0);
xhat = eta_l1(B*y,lam);
b=(1/M)*sum((abs(B*y)>lam).*(1-lam./(2*abs(B*y)+eps)),1);
c=(1/M)*sum((abs(B*y)>lam).*(0.5*lam.*sqrt(B*y).*(1./((conj(B*y)).^(3/2)+eps))),1);
% nmse(1,1)=10*log10(mean((sum(abs(xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))

for t=1:T-1
    switch t
        case 1
            theta=theta1(1);
            B=B_1;
        case 2
            theta=theta2(1);
            B=B_2;
        case 3
            theta=theta3(1);
            B=B_3;
        case 4
            theta=theta4(1);
            B=B_4;
        case 5
            theta=theta5(1);
            B=B_5;
        case 6
            theta=theta6(1);
            B=B_6;
        case 7
            theta=theta7(1);
            B=B_7;
    end
    v = y - A*xhat + bsxfun(@times,v,b) + bsxfun(@times,conj(v),c); % residual
    rhat = xhat + B*v; % denoiser input
    rvar=1/M*sum(abs(v).^2,1); % denoiser input err var
    lam=theta*sqrt(rvar);
    lam=max(abs(lam),0);
    xhat = eta_l1(rhat,lam); % estimate
	b=(1/M)*sum((abs(rhat)>lam).*(1-lam./(2*abs(rhat)+eps)),1);
    c=(1/M)*sum((abs(rhat)>lam).*(0.5*lam.*sqrt(rhat).*(1./((conj(rhat)).^(3/2)+eps))),1);
%     nmse(1,t+1)=10*log10(mean((sum(abs(xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
end
