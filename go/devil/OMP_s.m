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
        r = yl - A*xl;  %残差
        if (r'*r<th)
            break;
        end
        mf = A'*r;
        normA = diag(A'*A);
        [~,idx] = max(abs(mf./normA));  
        om = [om idx];
        Aom = A(:,om);  %找出投影列向量位置，并将多于列置0
        xl(om) = (Aom'*Aom)\(Aom'*yl);  %解出系数
    end
    x_est(:,l) = xl;  %构造系数矩阵
end

end

