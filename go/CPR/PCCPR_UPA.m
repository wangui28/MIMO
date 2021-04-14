function [h_est_list] = PCCPR_UPA(y, W, h0, sys, options)
%CRAF Summary of this function goes here
%   Detailed explanation goes here
h_est_list = zeros(sys.N,sys.T+1);

D1 = exp(-1i*2*pi*(0:1:(sys.N1-1))'*(0:1:(options.n1_dic-1))/options.n1_dic)/sqrt(sys.N1);
D2 = exp(-1i*2*pi*(0:1:(sys.N2-1))'*(0:1:(options.n2_dic-1))/options.n2_dic)/sqrt(sys.N2);
D = kron(D1, D2);  %������ɢ�������ֵ�
A = W'*D; %ϡ��������
h_sparse = D'*h0;  %ϡ���ŵ�ʸ��
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
    y_complex = y.*y_phase;  %��λ����

    residue_norm = norm(y_complex - W'*h_est);
    if (residue_norm < options.tol)
        for ii = iter:sys.T
            h_est_list(:,ii) = h_est;  %��ֹ����
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
%     end  %Ȩ��
%     h_sparse0 = mu*A'* diag(w)*(y_complex - A*h_sparse) + h_sparse;  %�ݶȷ�
%     [~, ii] = sort(h_sparse0, 'descend');
%     h_sparse = zeros(sys.N, 1);
%     h_sparse(ii(1:options.k)) = h_sparse0(ii(1:options.k));
%     h_est = D*h_sparse;
%     h_est_list(:,iter) = h_est;

    %% OMP 
    h_sparse = OMP(A, y_complex, options.k, 1e-4);  %OMP�㷨�ҵ����ϡ���ŵ���������
    h_est = D*h_sparse;  
    h_est_list(:,iter+1) = h_est;  %����ÿ�ν��

end

end

