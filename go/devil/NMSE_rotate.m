function [NMSE_list] = NMSE_rotate(h_est_list, h)
%NMSE_ROTATE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[~,T] = size(h_est_list);
NMSE_list = zeros(T,1);
for iter = 1:T
    h_est = squeeze(h_est_list(:,iter));  %������ǰ��
    alpha = (h_est'*h)/(h'*h);
    sol = alpha*h_est;
    NMSE_list(iter) = norm(sol - h, 2)^2;  %�����һ����С�������
end
end

