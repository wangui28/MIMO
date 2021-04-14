function [h, alpha, theta_t, theta_r] = Channel(Nr, Nt, L)
%CHANNEL 
%  信道建模
%     Nr：接收天线数
%     Nt：发送天线数
%     L：路径数

    lambda = 1; % 波长
    d = 0.5; % 天线间隔
    kappa = 10; % LoS/NLoS路径增益比值 
    L = random('Poisson',L)+1;
    
    %H = zeros(Nr, Nt);
    h = zeros(Nr*Nt, 1);
    alpha = zeros(L,1); % 路径增益 
    phi_t = zeros(L,1); % 离开角
    phi_r = zeros(L,1); % 到达角
    theta_t = zeros(L,1); % 离开方向矢量
    theta_r = zeros(L,1); % 到达方向矢量
    
    for l = 1:L
        phi_t(l) = 2*pi*rand(1);
        phi_r(l) = 2*pi*rand(1); % 角度范围[0, 2*pi)
        
        flag = 1;
        while (flag)
            for p = 1:l-1
                if (sqrt((phi_t(l)-phi_t(p))^2*Nt^2+(phi_r(l)-phi_r(p))^2*Nr^2)<1) % 控制角度间隔
                    phi_t(l) = 2*pi*rand(1);
                    phi_r(l) = 2*pi*rand(1); % 重新生成离开角、到达角
                    break;
                end    
            end
            flag = 0;
        end
    
        if (l==1)
            alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2)*sqrt(kappa); % 第1条LoS路径的归一化增益
        else
            alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2); % 第2—L条NLoS路径的归一化增益
        end
        
        theta_t(l) = 2*pi*sin(phi_t(l))*d/lambda;
        theta_r(l) = 2*pi*sin(phi_r(l))*d/lambda;
        
        at = 1/sqrt(Nt)*exp(1i*(0:1:Nt-1)'*theta_t(l)); % 离开方向矢量
        ar = 1/sqrt(Nr)*exp(1i*(0:1:Nr-1)'*theta_r(l)); % 到达方向矢量
        
        %H = H + alpha(l)*ar*at'; % 角度域信道矩阵
        %h = H(:);
        
        h = h + kron(at, ar)*alpha(l);
    end
    
    [alpha, alpha_sort] = sort(alpha, 'descend'); % 递减排序
    theta_t = theta_t(alpha_sort);
    theta_r = theta_r(alpha_sort); % 位置对齐
        
end

