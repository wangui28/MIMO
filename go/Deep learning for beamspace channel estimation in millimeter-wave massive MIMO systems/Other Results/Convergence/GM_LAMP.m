function[xhat,nmse]=GM_LAMP(y,A,trainedfile_name,T)
% y:measurement
% A:sensing matrix

[M,N]=size(A);

load(trainedfile_name);

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
        case 8
            theta=theta8;
            B=B_8;
        case 9
            theta=theta9;
            B=B_9;
        case 10
            theta=theta10;
            B=B_10;
        case 11
            theta=theta11;
            B=B_11;
    end
    v = y - A*xhat + bsxfun(@times,v,b) + bsxfun(@times,conj(v),c); % residual
    rhat = xhat + B*v; % denoiser input
    rvar=1/M*sum(abs(v).^2,1); % denoiser input err var
    [xhat,dxdr,dxdr1] = eta_gm(rhat,rvar,theta); % estimate
    b=(N/M)*dxdr;
    c=(N/M)*dxdr1;
%     nmse(1,t+1)=10*log10(mean((sum(abs(xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
end
