function [h_est_list] = PCCPR_SR_UPA(y, W, h0, sys, options)
%CRAF Summary of this function goes here
%   Detailed explanation goes here

%% initialization
% initial channel estimation (on grid)
h_est_list = zeros(sys.N,sys.T2+1); % N*(T2+1)
[h_est_list_grid] = PCCPR_UPA(y, W, h0, sys, options);
h_est = h_est_list_grid(:,sys.T+1);  %获得离散角度域信道估计
% h_est = h0;

% initial phase estimation
y_phase = y;
for i = 1:options.N_block
    W1 = W(:,(i-1)*sys.block_size+1:i*sys.block_size);
    y1 = y((i-1)*sys.block_size+1:i*sys.block_size);
    y1_phase = sign(y1'*(W1'*h_est));
    y_phase((i-1)*sys.block_size+1:i*sys.block_size) = y1_phase*ones(sys.block_size,1);
end
y_complex = y.*y_phase;  %相位补偿

% initial path coefficients estimation
[ h_est, om1, om2, alpha ] = NOMP_2D('1st', sys.N1, sys.N2, W, y_complex, options.Np, 1e-5, sys.T1);  %初始化路径参数
%disp(['init NMSE:',num2str(SE_rotate(h_est, h))]);  
h_est_list(:,1) = h_est;  %保存初始估计

%% iterative refinement
for iter = 1:sys.T2
    disp(['SR_iter=', num2str(iter)]);
    %% phase refinement based on the channel estimate of the previous iteration 
    y_phase = y;
    for i = 1:options.N_block
        W1 = W(:,(i-1)*sys.block_size+1:i*sys.block_size);
        y1 = y((i-1)*sys.block_size+1:i*sys.block_size);
        y1_phase = sign(y1'*(W1'*h_est));
        y_phase((i-1)*sys.block_size+1:i*sys.block_size) = y1_phase*ones(sys.block_size,1);
    end
    y_complex = y.*y_phase;  %相位补偿
    
    %% if residue is small enough, then break
    residue_norm = norm(y_complex - W'*h_est);
    if (residue_norm < options.tol)
        for ii = iter:(sys.T2+1)
            h_est_list(:,ii) = h_est;  
        end
        break;  
    end  %迭代终止判断
    disp(['residue_norm=', num2str(residue_norm)]);
    %% Gradient descend
    [h_est, om1, om2, alpha ] = Gradient_2D_once('1st', sys.N1, sys.N2, W, y_complex, om1, om2, alpha, h_est, 1e-6);%-1e-6/(20+(iter-1)*T_grad+loop)
    %disp(['iter:',num2str(iter),'   NMSE:',num2str(SE_rotate(h_est, h))]);
    h_est_list(:,iter+1) = h_est;  %梯度下降优化路径参数
    
end

end

