function [h0] = InitGuess_1(y, W, sys, options, var)
%% initial guess with the MM preprocessing function
D1 = exp(-1i*2*pi*(0:1:(sys.N1-1))'*(0:1:(options.n1_dic-1))/options.n1_dic)/sqrt(sys.N1);
D2 = exp(-1i*2*pi*(0:1:(sys.N2-1))'*(0:1:(options.n2_dic-1))/options.n2_dic)/sqrt(sys.N2);
D = kron(D1, D2);
A = W'*D;

corr = zeros(options.n1_dic*options.n2_dic, 1);
for block_idx = 1:options.N_block
    Ab = A((block_idx-1)*sys.block_size+1:block_idx*sys.block_size,:);
    Abn = Ab'*(inv(Ab*Ab'));
    corr = corr + abs(Abn*y((block_idx-1)*sys.block_size+1:block_idx*sys.block_size)).^2;  %֧����������ͳ����
end
[~, order] = sort(corr, 'descend');  %����������

if (options.init_supp == 1)
    Supp_opi = order(1);  %���Ԫ�ص�λ��
    z0 = zeros(options.n1_dic*options.n2_dic, 1);
    z0(Supp_opi) = 1;
    z0 = z0/norm(z0)*sqrt(y'*y/options.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
else
    supp = sort(order(1:options.init_supp));  %ȡ����ǰT��Ԫ�أ���ԭ������
    Supp_opi = sort(supp)';

% Supp_opi = find(x~=0);
% [~,aa] = sort(x(Supp_opi), 'descend');
% Supp_opi = Supp_opi(aa(1:k));

    ys = conj(y).*y;
    ys_mean = mean(ys);
    ys_normal = ys/ys_mean;
    delta = options.M/options.init_supp;
%   delta = options.M;
    ys_normal_plus = max(0, ys_normal);
    T = (ys_normal_plus-1)./(ys_normal_plus+sqrt(delta)-1);
    T = T*ys_mean;  %Ԥ������ֵ




    Aselectx = A(:,Supp_opi);  %����������
    D_mat = (1/options.M)*Aselectx' * diag(T) * (Aselectx);  %�����ʼ������
    [Vc, Dc] = eig(D_mat);  %����ֵ�ֽ⣬V������������D����ֵ�ԽǾ���
    [~, indeig]=sort(diag(abs(Dc)), 'descend');
    v0 = Vc(:, indeig(1));  %�ҳ��������ֵ�����Ӧ����������

    z0 = zeros(options.n1_dic*options.n2_dic, 1);
    z0(Supp_opi) = v0;
    z0 = z0/norm(z0)*sqrt(y'*y/options.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
end
% 
% x0 = zeros(n, 1);
% x0(Supp_opi) = x(Supp_opi);

end

