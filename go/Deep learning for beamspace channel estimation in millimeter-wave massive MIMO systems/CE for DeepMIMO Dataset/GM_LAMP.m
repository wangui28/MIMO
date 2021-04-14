function[xhat,nmse]=GM_LAMP(y,A,type,snr_range)
% y:measurement
% A:sensing matrix

[M,N]=size(A);
T=8; % the number of layer

if type==1 && snr_range==1
    load(['Trained_GM_LAMP_for_DeepMIMO_ULA_',num2str(N),num2str(M),'0to10dB.mat'])
end
if type==1 && snr_range==2
    load(['Trained_GM_LAMP_for_DeepMIMO_ULA_',num2str(N),num2str(M),'10to20dB.mat'])
end
if type==2 && snr_range==1
    load(['Trained_GM_LAMP_for_DeepMIMO_UPA_',num2str(N),num2str(M),'0to10dB.mat'])
end
if type==2 && snr_range==2
    load(['Trained_GM_LAMP_for_DeepMIMO_UPA_',num2str(N),num2str(M),'10to20dB.mat'])
end

v=y;  %Initialization of residual
theta=theta0;
B=B_0;
rvar = (1/M)*sum(abs(y).^2,1);
[xhat,dxdr,dxdr1] = eta_gm(B*y,rvar,theta);
b=(N/M)*dxdr;
c=(N/M)*dxdr1;
% nmse(1,1)=10*log10(mean((sum(abs(xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))

for t=1:T-1
    switch t
        case 1
            theta=theta1;
            B=B_1;
        case 2
            theta=theta2;
            B=B_2;
        case 3
            theta=theta3;
            B=B_3;
        case 4
            theta=theta4;
            B=B_4;
        case 5
            theta=theta5;
            B=B_5;
        case 6
            theta=theta6;
            B=B_6;
        case 7
            theta=theta7;
            B=B_7;
    end
    v = y - A*xhat + bsxfun(@times,v,b) + bsxfun(@times,conj(v),c); % residual
    rhat = xhat + B*v; % denoiser input
    rvar=1/M*sum(abs(v).^2,1); % denoiser input err var
    [xhat,dxdr,dxdr1] = eta_gm(rhat,rvar,theta); % estimate
    b=(N/M)*dxdr;
    c=(N/M)*dxdr1;
%     nmse(1,t+1)=10*log10(mean((sum(abs(xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
end
