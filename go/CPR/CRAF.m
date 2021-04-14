function [x_est, MSE] = CRAF( y, x, options, Amatrix, kB, B )
%CRAF

ymag  = sqrt(y);
mag_est = sqrt(sum(y)/options.m);


% based on squared quantities support recovery
rdata  = ((y)'*(abs(Amatrix).^2))'/options.m;    %支撑索引集的统计量
block_rdata =  zeros(round(options.n1/B),1);
for idx_b  = 1: round(options.n1/B)
    block_rdata(idx_b) = norm(rdata((B*idx_b-B+1):B*idx_b));
end
[~,sind] = sort(block_rdata, 'descend');   %构造索引集

block_supp = sind(1 : kB);
Sbhat = sort(block_supp);  
Supp_opi = [];
for supp_idx = 1:kB
    Supp_opi = [Supp_opi,(Sbhat(supp_idx)-1)*B+1:B*Sbhat(supp_idx)];  %取最大的前k个元素，按原序排列 
end
Supp_opi = Supp_opi';  
%% initilization
tau = 1; %M2 = 2;
%tweight = 3 * (y <= mag_est^2 ./ (1 + tau+1)) .* (-1)  + (1 ) .* (y > mag_est^2 ./ (1 + tau/2)); % changing the threshhold
tweight = 1.2 * (y <= mag_est^2 ./ (1 + 1)) .* (-1)  + (1 ) .* (y > mag_est^2 ./ (1 + 1)); % changing the threshhold权重更新
Aselectx = Amatrix(:,Supp_opi);
M = 1/options.m * Aselectx' * diag(tweight) * (Aselectx);  %初始化矩阵
[Vc, Dc] = eig(M);
Vr = real(Vc);
Dr = real(Dc);  %特征值分解
[~, indeig]=sort(diag(Dr), 'descend');
v0 = Vc(:, indeig(1));
 
z0 = zeros(options.n1, 1);
z0(Supp_opi) = v0;
x_est  =  z0/norm(z0)*sqrt(sum(y)/options.m);  % Apply scaling初始值

alpha = (x_est'*x)/(x'*x);
sol = alpha*x_est;
MSE = zeros(options.T,1);

for t = 1: options.T
    Az       = Amatrix * x_est; %A(z);
    %ratio    = abs(Az) ./ ymag;
    %yz       = ratio > 1 / (1 + Params.gamma_lb);
    ang      = options.cplx_flag * exp(1i * angle(Az)) + (1 - options.cplx_flag) * sign(Az);  %相位
    yz       = max(1 ./ ( 1 +  options.betaRAF* (ymag ./ abs(Az))), 0.1); % set the weight for reweighted amplitude flow迭代权重

    grad     = Amatrix' * (yz .* ymag .* ang - yz .* Az) / options.m;
    x_est        = x_est + options.muRAF * grad;  %梯度法搜寻新的估计值  
    
    block_z  = zeros(round(options.n1/B), 1);
    for idx_b = 1:round(options.n1/B)
        block_z(idx_b) =  norm(x_est((idx_b*B-B+1):idx_b*B), 2);
    end    
    [sz, sz_ind]  = sort(block_z, 'descend');
    S0 =  sz_ind(1:kB);
    block_idx = sort(S0);
    supp_z = [];
    for supp_idx = 1:kB
        supp_z = [supp_z, (block_idx(supp_idx)-1)*B+1:(block_idx(supp_idx)*B)];
    end
    non_supp = setdiff((1:options.n1), supp_z);
    x_est(non_supp) = 0;  %保留最大元素并将其他元素置0
    
    alpha = (x_est'*x)/(x'*x);
    sol = alpha*x_est;
    MSE(t) = norm(sol - x, 2)^2;  %计算均方误差
    if (MSE(t) < options.tol)
        MSE(t+1:end) = MSE(t+1:end)+MSE(t);
        break;
    end
end
end

