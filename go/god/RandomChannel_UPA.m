function [ h, psi, theta, alpha ] = RandomChannel_UPA(system)
%RANDOMCHANNEL 
%   �ŵ���ģ

L = random('Poisson',system.L)+1;
% L = sys.L;
psi = zeros(L,1);  %���Ǿ���
theta = zeros(L,1);  %��λ�Ǿ���
alpha = zeros(L,1);  %·���������
h = zeros(system.N,1);  %�����ŵ�����

for l = 1:L
    psi(l) = pi*rand(1)-pi/2;
    theta(l) = pi*rand(1)-pi/2;
    
    %·������
    flag = 1;
    while (flag)
        for pp = 1:l-1
            if (sqrt((psi(l)-psi(pp))^2*system.N1*system.N1+(theta(l)-theta(pp))^2*system.N2*system.N2)<1)  %���ƽǶȼ��
                psi(l) = pi*rand(1)-pi/2;
                theta(l) = pi*rand(1)-pi/2;  %��������·������
                break;
            end    
        end
        flag = 0;
    end
    
    %·������
    if (l==1)
        alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2)*sqrt(system.kappa);  %��1��LoS·���Ĺ�һ������
    else
        alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2);  %��2 -- L��LoS·���Ĺ�һ������
    end
    operator = 1/sqrt(L-1+system.kappa);  %����·����������
    
    
    h = h + 1/sqrt(system.N)*operator*alpha(l)*kron(exp(1i*pi*sin(psi(l))*(0:1:(system.N1-1)))', exp(1i*pi*cos(psi(l))*sin(theta(l))*(0:1:(system.N2-1)))');
    %�����ŵ�ģ��N*1
end

% h = h/norm(h);
% normh = norm(h,2);
% h = sqrt(N)*h/normh;
% alpha = sqrt(N)*alpha/normh;
alpha = alpha.*operator;  %�ŵ�����
% h = h/sqrt(sys.N);
[alpha, alpha_sort] = sort(alpha, 'descend');  %�ݼ�����
psi = psi(alpha_sort);
theta = theta(alpha_sort);  %λ�ö���
end

