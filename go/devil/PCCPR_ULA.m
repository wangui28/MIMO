function [h_est_list] = PCCPR_ULA(y, C, h0, sys)
%PCCPR_ULA 
%   

h_est_list = zeros(sys.N,sys.T+1);

D1 = exp(-1i*2*pi*(0:1:(sys.Nt-1))'*(0:1:(sys.Nt-1))/(sys.Nt-1))/sqrt(sys.Nt); %�Ƕ����뿪����ʸ��
D2 = exp(-1i*2*pi*(0:1:(sys.Nr-1))'*(0:1:(sys.Nr-1))/(sys.Nr-1))/sqrt(sys.Nr); %�Ƕ��򵽴﷽��ʸ��
D = kron(conj(D1), D2); %�Ƕ����ֵ�
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
    y_complex = y.*y_phase;  %��λ����
    
    residue_norm = norm(y_complex - C*h_est);
    if (residue_norm < sys.tol)
        for j = iter:sys.T
            h_est_list(:,j) = h_est;
        end
        break;
    end %��ֹ����Ӳ��ֵ�ж�
    
%    %% reweighted HT
%     mu = 0.2;
%     tao_w = 0.1;
%     beta = 0.1;
%     weight = zeros(sys.M,1);
%     for k = 1:sys.M
%         weight(k) = max(tao_w, abs(A(k,:)*h_sparse)/(abs(A(k,:)*h_sparse)+beta*abs(y(k))));
%     end  %Ȩ�ظ���
%     h_sparse0 = mu*A'* diag(weight)*(y - A*h_sparse) + h_sparse;  %�ݶ��ؼ�Ȩ
%     [~,j] = sort(h_sparse0, 'descend');
%     h_sparse = zeros(sys.N, 1);
%     h_sparse(j(1:sys.k)) = h_sparse0(j(1:sys.k)); %ѡȡk�����Ԫ��
%     h_est = D*h_sparse;
%     h_est_list(:,iter) = h_est;

    %% OMP 
    h_sparse = OMP_s(A, y_complex, sys.k, 1e-4);  %OMP�㷨�ҵ���ѽǶ����ŵ���������
    h_est = D*h_sparse;  
    h_est_list(:,iter+1) = h_est;  %����ÿ�ν��

end
end
