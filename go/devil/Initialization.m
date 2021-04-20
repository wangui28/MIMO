function [h0] = Initialization(y, C, sys, alg)
%initialization
%   ��ʼ���ŵ�����ֵh0
D1 = exp(-1i*2*pi*(0:1:(sys.Nt-1))'*(0:1:(alg.Nt-1))/(alg.Nt-1))/sqrt(sys.Nt); %�Ƕ����뿪����ʸ��
D2 = exp(-1i*2*pi*(0:1:(sys.Nr-1))'*(0:1:(alg.Nr-1))/(alg.Nr-1))/sqrt(sys.Nr); %�Ƕ��򵽴﷽��ʸ��
D = kron(D1, D2); %�Ƕ����ֵ�
A = C'*D;

corr = zeros(alg.Nt*alg.Nr, 1);

for block_idx = 1:alg.block_num
    Ab = A((block_idx-1)*alg.block_size+1:block_idx*alg.block_size, :);
    Abn = Ab'*(inv(Ab*Ab'));
    corr = corr + abs(Abn*y((block_idx-1)*alg.block_size+1:block_idx*alg.block_size)).^2;  %֧����������ͳ����
end
[~, order] = sort(corr, 'descend');  %����������

if (alg.init_supp == 1)  %������1��
    Supp_opi = order(1);  %���Ԫ�ص�λ��
    z0 = zeros(alg.Nt*alg.Nr, 1);
    z0(Supp_opi) = 1;
    z0 = z0/norm(z0)*sqrt(y'*y/alg.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
else
    supp = sort(order(1:alg.init_supp));  %ȡ����ǰinit_supp��Ԫ�أ���ԭ������
    Supp_opi = sort(supp)';

    ys = conj(y).*y;
    ys_mean = mean(ys);
    T = zeros(alg.M,1);
    I_plus = find(ys>(ys_mean/2));
    I_minus = find(ys<=(ys_mean/2));  %����������
    T(I_plus) = 1/numel(I_plus);
    T(I_minus) = -1/numel(I_minus);  %Ԥ������ֵ

    Aselectx = A(:,Supp_opi);  %����������
    D_mat = (1/alg.M)*Aselectx' * diag(T) * Aselectx;  %�����ʼ������
    [Vc, Dc] = eig(D_mat);  %����ֵ�ֽ⣬V������������D����ֵ�ԽǾ���
    [~, indeig]=sort(diag(abs(Dc)), 'descend');
    v0 = Vc(:, indeig(1));  %�ҳ��������ֵ�����Ӧ����������

    z0 = zeros(alg.Nt*alg.Nr, 1);
    z0(Supp_opi) = v0;
    z0 = z0/norm(z0)*sqrt(y'*y/alg.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
end

end
