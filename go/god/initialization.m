function [h0] = initialization(y,W,system,T_max)
%INITIALIZATION
%   信道初值
D1 = exp(-1i*2*pi*(0:1:(system.N1-1))'*(0:1:(system.N1-1))/system.N1)/sqrt(system.N1);  %尺寸N1*N1的DFT
D2 = exp(-1i*2*pi*(0:1:(system.N2-1))'*(0:1:(system.N2-1))/system.N2)/sqrt(system.N2);  %尺寸N2*N2的DFT
D = kron(D1, D2);  %构造离散域网格 N*N
A = W'*D;  %离散组合矩阵 M*N

corr = zeros(system.N1*system.N2, 1);

for block_idx = 1:system.N_block
    Ab = A((block_idx-1)* system.N_RF+1:block_idx* system.N_RF,:);
    Abn = Ab'*(inv(Ab*Ab'));
    corr = corr + abs(Abn*y((block_idx-1)* system.N_RF+1:block_idx* system.N_RF)).^2;  %支撑索引集的统计量
end
[~, order] = sort(corr, 'descend');  %构造索引集

if (T_max == 1)  %仅迭代1次
    Supp_opi = order(1);  %最大元素的位置
    z0 = zeros(system.N1*system.N1, 1);
    z0(Supp_opi) = 1;
    z0 = z0/norm(z0)*sqrt(y'*y/system.M);
    h0 = D*z0;  %生成初始角度域信道向量
else
    supp = sort(order(1:T_max));  %取最大的前T个元素，按原序排列
    Supp_opi = sort(supp)';
    
    ys = conj(y).*y;
    ys_mean = mean(ys);
    Tao = zeros(system.M,1);
    I_plus = find(ys>(ys_mean/2));
    I_minus = find(ys<=(ys_mean/2));  %索引集分组
    Tao(I_plus) = 1/numel(I_plus);
    Tao(I_minus) = -1/numel(I_minus);  %预处理函数值
    
    Aselect = A(:,Supp_opi);  %索引列向量
    D_mat = (1/system.M)*Aselect' * diag(Tao) * (Aselect);  %构造初始化矩阵
    [Vc, Dc] = eig(D_mat);  %特征值分解，V特征向量矩阵、D特征值对角矩阵
    [~, indeig]=sort(diag(abs(Dc)), 'descend');
    v0 = Vc(:, indeig(1));  %找出最大特征值及其对应的特征向量

    z0 = zeros(system.N1*system.N2, 1);
    z0(Supp_opi) = v0;
    z0 = z0/norm(z0)*sqrt(y'*y/system.M);
    h0 = D*z0;  %生成初始角度域信道向量
end

end

