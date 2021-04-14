function [h0] = initialization(y,W,system,T_max)
%INITIALIZATION
%   �ŵ���ֵ
D1 = exp(-1i*2*pi*(0:1:(system.N1-1))'*(0:1:(system.N1-1))/system.N1)/sqrt(system.N1);  %�ߴ�N1*N1��DFT
D2 = exp(-1i*2*pi*(0:1:(system.N2-1))'*(0:1:(system.N2-1))/system.N2)/sqrt(system.N2);  %�ߴ�N2*N2��DFT
D = kron(D1, D2);  %������ɢ������ N*N
A = W'*D;  %��ɢ��Ͼ��� M*N

corr = zeros(system.N1*system.N2, 1);

for block_idx = 1:system.N_block
    Ab = A((block_idx-1)* system.N_RF+1:block_idx* system.N_RF,:);
    Abn = Ab'*(inv(Ab*Ab'));
    corr = corr + abs(Abn*y((block_idx-1)* system.N_RF+1:block_idx* system.N_RF)).^2;  %֧����������ͳ����
end
[~, order] = sort(corr, 'descend');  %����������

if (T_max == 1)  %������1��
    Supp_opi = order(1);  %���Ԫ�ص�λ��
    z0 = zeros(system.N1*system.N1, 1);
    z0(Supp_opi) = 1;
    z0 = z0/norm(z0)*sqrt(y'*y/system.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
else
    supp = sort(order(1:T_max));  %ȡ����ǰT��Ԫ�أ���ԭ������
    Supp_opi = sort(supp)';
    
    ys = conj(y).*y;
    ys_mean = mean(ys);
    Tao = zeros(system.M,1);
    I_plus = find(ys>(ys_mean/2));
    I_minus = find(ys<=(ys_mean/2));  %����������
    Tao(I_plus) = 1/numel(I_plus);
    Tao(I_minus) = -1/numel(I_minus);  %Ԥ������ֵ
    
    Aselect = A(:,Supp_opi);  %����������
    D_mat = (1/system.M)*Aselect' * diag(Tao) * (Aselect);  %�����ʼ������
    [Vc, Dc] = eig(D_mat);  %����ֵ�ֽ⣬V������������D����ֵ�ԽǾ���
    [~, indeig]=sort(diag(abs(Dc)), 'descend');
    v0 = Vc(:, indeig(1));  %�ҳ��������ֵ�����Ӧ����������

    z0 = zeros(system.N1*system.N2, 1);
    z0(Supp_opi) = v0;
    z0 = z0/norm(z0)*sqrt(y'*y/system.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
end

end

