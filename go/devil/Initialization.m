function [h0] = Initialization(y, C, sys, alg)
%initialization
%   初始化信道估计值h0
D1 = exp(-1i*2*pi*(0:1:(sys.Nt-1))'*(0:1:(alg.Nt-1))/(alg.Nt-1))/sqrt(sys.Nt); %角度域离开方向矢量
D2 = exp(-1i*2*pi*(0:1:(sys.Nr-1))'*(0:1:(alg.Nr-1))/(alg.Nr-1))/sqrt(sys.Nr); %角度域到达方向矢量
D = kron(D1, D2); %角度域字典
A = C'*D;

corr = zeros(alg.Nt*alg.Nr, 1);

for block_idx = 1:alg.block_num
    Ab = A((block_idx-1)*alg.block_size+1:block_idx*alg.block_size, :);
    Abn = Ab'*(inv(Ab*Ab'));
    corr = corr + abs(Abn*y((block_idx-1)*alg.block_size+1:block_idx*alg.block_size)).^2;  %支撑索引集的统计量
end
[~, order] = sort(corr, 'descend');  %构造索引集

if (alg.init_supp == 1)  %仅迭代1次
    Supp_opi = order(1);  %最大元素的位置
    z0 = zeros(alg.Nt*alg.Nr, 1);
    z0(Supp_opi) = 1;
    z0 = z0/norm(z0)*sqrt(y'*y/alg.M);
    h0 = D*z0;  %生成初始角度域信道向量
else
    supp = sort(order(1:alg.init_supp));  %取最大的前init_supp个元素，按原序排列
    Supp_opi = sort(supp)';

    ys = conj(y).*y;
    ys_mean = mean(ys);
    T = zeros(alg.M,1);
    I_plus = find(ys>(ys_mean/2));
    I_minus = find(ys<=(ys_mean/2));  %索引集分组
    T(I_plus) = 1/numel(I_plus);
    T(I_minus) = -1/numel(I_minus);  %预处理函数值

    Aselectx = A(:,Supp_opi);  %索引列向量
    D_mat = (1/alg.M)*Aselectx' * diag(T) * Aselectx;  %构造初始化矩阵
    [Vc, Dc] = eig(D_mat);  %特征值分解，V特征向量矩阵、D特征值对角矩阵
    [~, indeig]=sort(diag(abs(Dc)), 'descend');
    v0 = Vc(:, indeig(1));  %找出最大特征值及其对应的特征向量

    z0 = zeros(alg.Nt*alg.Nr, 1);
    z0(Supp_opi) = v0;
    z0 = z0/norm(z0)*sqrt(y'*y/alg.M);
    h0 = D*z0;  %生成初始角度域信道向量
end

end
