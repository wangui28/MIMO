function [ h, psi, theta, alpha ] = RandomChannel_UPA(system)
%RANDOMCHANNEL 
%   信道建模

L = random('Poisson',system.L)+1;
% L = sys.L;
psi = zeros(L,1);  %仰角矩阵
theta = zeros(L,1);  %方位角矩阵
alpha = zeros(L,1);  %路径增益矩阵
h = zeros(system.N,1);  %生成信道矩阵

for l = 1:L
    psi(l) = pi*rand(1)-pi/2;
    theta(l) = pi*rand(1)-pi/2;
    
    %路径方向
    flag = 1;
    while (flag)
        for pp = 1:l-1
            if (sqrt((psi(l)-psi(pp))^2*system.N1*system.N1+(theta(l)-theta(pp))^2*system.N2*system.N2)<1)  %控制角度间隔
                psi(l) = pi*rand(1)-pi/2;
                theta(l) = pi*rand(1)-pi/2;  %重新生成路径方向
                break;
            end    
        end
        flag = 0;
    end
    
    %路径增益
    if (l==1)
        alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2)*sqrt(system.kappa);  %第1条LoS路径的归一化增益
    else
        alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2);  %第2 -- L条LoS路径的归一化增益
    end
    operator = 1/sqrt(L-1+system.kappa);  %定义路径增益算子
    
    
    h = h + 1/sqrt(system.N)*operator*alpha(l)*kron(exp(1i*pi*sin(psi(l))*(0:1:(system.N1-1)))', exp(1i*pi*cos(psi(l))*sin(theta(l))*(0:1:(system.N2-1)))');
    %构造信道模型N*1
end

% h = h/norm(h);
% normh = norm(h,2);
% h = sqrt(N)*h/normh;
% alpha = sqrt(N)*alpha/normh;
alpha = alpha.*operator;  %信道增益
% h = h/sqrt(sys.N);
[alpha, alpha_sort] = sort(alpha, 'descend');  %递减排序
psi = psi(alpha_sort);
theta = theta(alpha_sort);  %位置对齐
end

