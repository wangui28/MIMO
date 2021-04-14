function[xhat]=AMP(y,A)

% y: measurement
% A: sensing matrix

[M,N]=size(A);
S=size(y,2);  
T=10; % number of iteration
alf=1.1402; % shrinkage parameter
eta_l1=@(r,lam) exp(1i*angle(r)).*max(bsxfun(@minus,abs(r),lam),0); %soft-thresholding shrinkage function

Bmf=A'; % matched filter
xhat=zeros(N,S); % Initialization of signal estimation
v=y; % Initialization of residual
% nmse=zeros(1,T);
for t=0:T-1
    rhat=xhat+Bmf*v;
    rvar=sum(abs(v).^2,1)/M;
    lam=alf*sqrt(rvar);
    xhat=eta_l1(rhat,lam);
	g=(1/M)*sum((xhat~=0).*(1-lam./(2*abs(rhat)+eps)),1);
    c=(1/M)*sum((xhat~=0).*(0.5*lam.*sqrt(rhat).*(1./((conj(rhat)).^(3/2)+eps))),1);
	v=y-A*xhat+bsxfun(@times,v,g)+bsxfun(@times,conj(v),c);
%     nmse(1,t+1)=10*log10(mean((sum(abs(xhat-x).^2,1)))./mean(sum(abs(x).^2,1)));
end
