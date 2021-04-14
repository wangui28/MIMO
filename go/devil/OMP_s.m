function [ x_est ] = OMP_s( A, y, T, th)
%OMP Summary of this function goes here
% Detailed explanation goes here
[~, N] = size(A);
[~, L] = size(y); % L=1
x_est = zeros(N,L);
for l = 1:L
    yl = y(:,l);  %y
    xl = zeros(N,1);
    om = [];
    for t = 1:T
        r = yl - A*xl;  %�в�
        if (r'*r<th)
            break;
        end
        mf = A'*r;
        normA = diag(A'*A);
        [~,idx] = max(abs(mf./normA));  
        om = [om idx];
        Aom = A(:,om);  %�ҳ�ͶӰ������λ�ã�������������0
        xl(om) = (Aom'*Aom)\(Aom'*yl);  %���ϵ��
    end
    x_est(:,l) = xl;  %����ϵ������
end

end

