function [h_est_list] = PCCPR_UPA(y, W, h0, sys, options)
%CRAF Summary of this function goes here
%   Detailed explanation goes here
h_est_list = zeros(sys.N,sys.T+1);

D1 = exp(-1i*2*pi*(0:1:(sys.N1-1))'*(0:1:(options.n1_dic-1))/options.n1_dic)/sqrt(sys.N1);
D2 = exp(-1i*2*pi*(0:1:(sys.N2-1))'*(0:1:(options.n2_dic-1))/options.n2_dic)/sqrt(sys.N2);
D = kron(D1, D2);  %构造离散域网络字典
A = W'*D; %稀疏编码矩阵
h_sparse = D'*h0;  %稀疏信道矢量
h_est = h0;
h_est_list(:,1) = h0;
for iter = 1:sys.T
    y_phase = y;
    for i = 1:options.N_block
        W1 = W(:,(i-1)*sys.block_size+1:i*sys.block_size);
        y1 = y((i-1)*sys.block_size+1:i*sys.block_size);
        y1_phase = sign(y1'*(W1'*h_est));
        y_phase((i-1)*sys.block_size+1:i*sys.block_size) = y1_phase*ones(sys.block_size,1);
    end
    y_complex = y.*y_phase;  %相位补偿

    residue_norm = norm(y_complex - W'*h_est);
    if (residue_norm < options.tol)
        for ii = iter:sys.T
            h_est_list(:,ii) = h_est;  %终止迭代
        end
        break;
    end
    
    
%    %% reweighted HT
%     mu = 0.2;
%     taow = 0.1;
%     beta = 0.1;
%     w = zeros(sys.M,1);
%     for ww = 1:sys.M
%         w(ww) = max(taow, abs(A(ww,:)*h_sparse)/(abs(A(ww,:)*h_sparse)+beta*abs(y(ww))));
%     end  %权重
%     h_sparse0 = mu*A'* diag(w)*(y_complex - A*h_sparse) + h_sparse;  %梯度法
%     [~, ii] = sort(h_sparse0, 'descend');
%     h_sparse = zeros(sys.N, 1);
%     h_sparse(ii(1:options.k)) = h_sparse0(ii(1:options.k));
%     h_est = D*h_sparse;
%     h_est_list(:,iter) = h_est;

    %% OMP 
    h_sparse = OMP(A, y_complex, options.k, 1e-4);  %OMP算法找到最佳稀疏信道向量估计
    h_est = D*h_sparse;  
    h_est_list(:,iter+1) = h_est;  %保留每次结果

end

end

