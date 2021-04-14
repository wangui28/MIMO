function [h_est_list] = PCCPR_ULA(y, C, h0, sys)
%PCCPR_ULA 
%   

h_est_list = zeros(sys.N,sys.T+1);

D1 = exp(-1i*2*pi*(0:1:(sys.Nt-1))'*(0:1:(sys.Nt-1))/(sys.Nt-1))/sqrt(sys.Nt); %角度域离开方向矢量
D2 = exp(-1i*2*pi*(0:1:(sys.Nr-1))'*(0:1:(sys.Nr-1))/(sys.Nr-1))/sqrt(sys.Nr); %角度域到达方向矢量
D = kron(conj(D1), D2); %角度域字典
A = C*D;

h_est = h0;
h_est_list(:,1) = h0;

for iter = 1:sys.T
    
	y_phase = y;
for i = 1:sys.block
    C1 = C((i-1)*sys.size+1:i*sys.size, :);
    y1 = y((i-1)*sys.size+1:i*sys.size);
    y1_phase = sign(y1'*(C1*h_est));
    y_phase((i-1)*sys.size+1:i*sys.size) = y1_phase*ones(sys.size,1);
end
    y_complex = y.*y_phase;  %相位补偿
    
    residue_norm = norm(y_complex - C*h_est);
    if (residue_norm < sys.tol)
        for j = iter:sys.T
            h_est_list(:,j) = h_est;
        end
        break;
    end %终止迭代硬阈值判断
    
%    %% reweighted HT
%     mu = 0.2;
%     tao_w = 0.1;
%     beta = 0.1;
%     weight = zeros(sys.M,1);
%     for k = 1:sys.M
%         weight(k) = max(tao_w, abs(A(k,:)*h_sparse)/(abs(A(k,:)*h_sparse)+beta*abs(y(k))));
%     end  %权重更新
%     h_sparse0 = mu*A'* diag(weight)*(y - A*h_sparse) + h_sparse;  %梯度重加权
%     [~,j] = sort(h_sparse0, 'descend');
%     h_sparse = zeros(sys.N, 1);
%     h_sparse(j(1:sys.k)) = h_sparse0(j(1:sys.k)); %选取k个最大元素
%     h_est = D*h_sparse;
%     h_est_list(:,iter) = h_est;

    %% OMP 
    h_sparse = OMP_s(A, y_complex, sys.k, 1e-4);  %OMP算法找到最佳角度域信道向量估计
    h_est = D*h_sparse;  
    h_est_list(:,iter+1) = h_est;  %保留每次结果

end
end

