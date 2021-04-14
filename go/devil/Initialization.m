function [h0] = Initialization(y, C, sys)
%INITIALIZATION
%   ��ʼ���ŵ�����ֵh0
D1 = exp(1i*2*pi*(0:1:(sys.Nt-1))'*(0:1:(sys.Nt-1))/(sys.Nt-1))/sqrt(sys.Nt); %�Ƕ����뿪����ʸ��
D2 = exp(1i*2*pi*(0:1:(sys.Nr-1))'*(0:1:(sys.Nr-1))/(sys.Nr-1))/sqrt(sys.Nr); %�Ƕ��򵽴﷽��ʸ��
D = kron(conj(D1), D2); %�Ƕ����ֵ�
A = C*D;

corr = zeros(sys.Nt*sys.Nr, 1);

for block_idx = 1:sys.block
    Ab = A((block_idx-1)*sys.size+1:block_idx*sys.size,:);
    Abn = Ab'*(inv(Ab*Ab'));
    corr = corr + abs(Abn*y((block_idx-1)*sys.size+1:block_idx*sys.size)).^2;  %֧����������ͳ����
end
[~, order] = sort(corr, 'descend');  %����������

if (sys.init_supp == 1)  %������1��
    Supp_opi = order(1);  %���Ԫ�ص�λ��
    z0 = zeros(sys.Nt*sys.Nr, 1);
    z0(Supp_opi) = 1;
    z0 = z0/norm(z0)*sqrt(y'*y/sys.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
else
    supp = sort(order(1:sys.init_supp));  %ȡ����ǰinit_supp��Ԫ�أ���ԭ������
    Supp_opi = sort(supp)';

% Supp_opi = find(x~=0);
% [~,aa] = sort(x(Supp_opi), 'descend');
% Supp_opi = Supp_opi(aa(1:k));

    ys = conj(y).*y;
    ys_mean = mean(ys);
    T = zeros(sys.M,1);
    I_plus = find(ys>(ys_mean/2));
    I_minus = find(ys<=(ys_mean/2));  %����������
    T(I_plus) = 1/numel(I_plus);
    T(I_minus) = -1/numel(I_minus);  %Ԥ��������ֵ
    
%     ys = conj(y).*y;
%     ys_mean = mean(ys);
%     ys_normal = ys/ys_mean;
%     delta = sys.M/sys.init_supp;
%     delta = options.M;
%     ys_normal_plus = max(0, ys_normal);
%     T = (ys_normal_plus-1)./(ys_normal_plus+sqrt(delta)-1);
%     T = T*ys_mean;

%     y = conj(y).*y;
%     sigma = sqrt(var*2);
%     eta_0 = 1/2 - exp(sigma.^2/2).*normcdf(-sigma);
%     mu_0 = 1/2 + (sigma.^2 - 1).*exp(sigma.^2/2).*normcdf(-sigma) - sigma/sqrt(2*pi);
%     T = (1 - (y-sigma.^2+sigma.*normpdf(y./sigma-sigma)./normcdf(y./sigma-sigma)).^(-1));
%     for nn = 1:sys.M
%         if (y(nn)<0)
%             T(nn) = 1-eta_0/mu_0;
%         end
%     end
    
    Aselectx = A(:,Supp_opi);  %����������
    D_mat = (1/sys.M)*Aselectx' * diag(T) * Aselectx;  %�����ʼ������
    [Vc, Dc] = eig(D_mat);  %����ֵ�ֽ⣬V������������D����ֵ�ԽǾ���
    [~, indeig]=sort(diag(abs(Dc)), 'descend');
    v0 = Vc(:, indeig(1));  %�ҳ��������ֵ�����Ӧ����������

    z0 = zeros(sys.Nt*sys.Nr, 1);
    z0(Supp_opi) = v0;
    z0 = z0/norm(z0)*sqrt(y'*y/sys.M);
    h0 = D*z0;  %���ɳ�ʼ�Ƕ����ŵ�����
end
% 
% x0 = zeros(n, 1);
% x0(Supp_opi) = x(Supp_opi);
end