function [h_est_list] = PCCPR_SR_ULA(y, C, h0, sys)
%CRAF Summary of this function goes here
%   Detailed explanation goes here

%% initialization
% initial channel estimation (on grid)
h_est_list = zeros(sys.N,sys.T2+1); % N*(T2+1)
[h_est_list_grid] = PCCPR_ULA(y, C, h0, sys);
h_est = h_est_list_grid(:,sys.T+1);  %获得离散角度域信道估计
% y = C*h_est;

y_phase = y;
for i = 1:sys.block
    C1 = C((i-1)*sys.size+1:i*sys.size, :);
    y1 = y((i-1)*sys.size+1:i*sys.size);
    y1_phase = sign(y1'*(C1*h_est));
    y_phase((i-1)*sys.size+1:i*sys.size) = y1_phase*ones(sys.size,1);
end
y_complex = y.*y_phase;  %相位补偿

% initial path coefficients estimation
[ h_est, theta_t, thera_r, alpha] = NOMP('1st', sys.Nt, sys.Nr, C, y_complex,  sys.L, 1e-5, sys.T1);  %初始化路径参数
%disp(['init NMSE:',num2str(SE_rotate(h_est, h))]);  
h_est_list(:,1) = h_est;  %保存初始估计

%% iterative refinement
for iter = 1:sys.T2
    disp(['SR_iter=', num2str(iter)]);
    y_phase = y;
for i = 1:sys.block
    C1 = C((i-1)*sys.size+1:i*sys.size, :);
    y1 = y((i-1)*sys.size+1:i*sys.size);
    y1_phase = sign(y1'*(C1*h_est));
    y_phase((i-1)*sys.size+1:i*sys.size) = y1_phase*ones(sys.size,1);
end
    y_complex = y.*y_phase;  %相位补偿

    %% if residue is small enough, then break
    residue_norm = norm(y_complex - C*h_est);
    if (residue_norm < sys.tol)
        for ii = iter:(sys.T2+1)
            h_est_list(:,ii) = h_est;  
        end
        break;  
    end  %终止迭代硬阈值判断
    disp(['residue_norm=', num2str(residue_norm)]);
    
    %% Gradient descend
    [h_est, theta_t, thera_r, alpha] = Grad_desc('1st', sys.Nt, sys.Nr, C, y_complex, theta_t, thera_r, alpha, h_est, 1e-6);%-1e-6/(20+(iter-1)*T_grad+loop)
    %disp(['iter:',num2str(iter),'   NMSE:',num2str(SE_rotate(h_est, h))]);
    h_est_list(:,iter+1) = h_est;  %梯度下降优化路径参数
   
end

end

