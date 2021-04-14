function [ SE_list ] = SE_rotate(h_est_list, h)
%NMSE_ROTATE Summary of this function goes here
%   Detailed explanation goes here
[~,T] = size(h_est_list);
SE_list = zeros(T,1);
for iter = 1:T
    h_est = squeeze(h_est_list(:,iter));  %保留当前列
    alpha = (h_est'*h)/(h'*h);
    sol = alpha*h_est;
    SE_list(iter) = norm(sol - h, 2)^2;  %计算最小均方误差
end
end

