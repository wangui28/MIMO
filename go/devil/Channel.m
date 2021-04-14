function [h, alpha, theta_t, theta_r] = Channel(Nr, Nt, L)
%CHANNEL 
%  �ŵ���ģ
%     Nr������������
%     Nt������������
%     L��·����

    lambda = 1; % ����
    d = 0.5; % ���߼��
    kappa = 10; % LoS/NLoS·�������ֵ 
    L = random('Poisson',L)+1;
    
    %H = zeros(Nr, Nt);
    h = zeros(Nr*Nt, 1);
    alpha = zeros(L,1); % ·������ 
    phi_t = zeros(L,1); % �뿪��
    phi_r = zeros(L,1); % �����
    theta_t = zeros(L,1); % �뿪����ʸ��
    theta_r = zeros(L,1); % ���﷽��ʸ��
    
    for l = 1:L
        phi_t(l) = 2*pi*rand(1);
        phi_r(l) = 2*pi*rand(1); % �Ƕȷ�Χ[0, 2*pi)
        
        flag = 1;
        while (flag)
            for p = 1:l-1
                if (sqrt((phi_t(l)-phi_t(p))^2*Nt^2+(phi_r(l)-phi_r(p))^2*Nr^2)<1) % ���ƽǶȼ��
                    phi_t(l) = 2*pi*rand(1);
                    phi_r(l) = 2*pi*rand(1); % ���������뿪�ǡ������
                    break;
                end    
            end
            flag = 0;
        end
    
        if (l==1)
            alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2)*sqrt(kappa); % ��1��LoS·���Ĺ�һ������
        else
            alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2); % ��2��L��NLoS·���Ĺ�һ������
        end
        
        theta_t(l) = 2*pi*sin(phi_t(l))*d/lambda;
        theta_r(l) = 2*pi*sin(phi_r(l))*d/lambda;
        
        at = 1/sqrt(Nt)*exp(1i*(0:1:Nt-1)'*theta_t(l)); % �뿪����ʸ��
        ar = 1/sqrt(Nr)*exp(1i*(0:1:Nr-1)'*theta_r(l)); % ���﷽��ʸ��
        
        %H = H + alpha(l)*ar*at'; % �Ƕ����ŵ�����
        %h = H(:);
        
        h = h + kron(at, ar)*alpha(l);
    end
    
    [alpha, alpha_sort] = sort(alpha, 'descend'); % �ݼ�����
    theta_t = theta_t(alpha_sort);
    theta_r = theta_r(alpha_sort); % λ�ö���
        
end

