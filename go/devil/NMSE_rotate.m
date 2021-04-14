function [NMSE_list] = NMSE_rotate(h_est_list, h)
%NMSE_ROTATE 此处显示有关此函数的摘要
%   此处显示详细说明
[~,T] = size(h_est_list);
NMSE_list = zeros(T,1);
for iter = 1:T
    h_est = squeeze(h_est_list(:,iter));  %保留当前列
    alpha = (h_est'*h)/(h'*h);
    sol = alpha*h_est;
    NMSE_list(iter) = norm(sol - h, 2)^2;  %计算归一化最小均方误差
end
end

