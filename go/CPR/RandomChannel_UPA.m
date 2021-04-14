function [ h, psi, theta, alpha ] = RandomChannel_UPA(sys)
%RANDOMCHANNEL Summary of this function goes here
%   Detailed explanation goes here
L = random('Poisson',sys.L)+1;
% L = sys.L;
psi = zeros(L,1);  %�������Ǿ���
theta = zeros(L,1);  %���ɷ�λ�Ǿ���
alpha = zeros(L,1);  %����·���������
h = zeros(sys.N,1);  %�����ŵ�����

for l = 1:L
    psi(l) = pi*rand(1)-pi/2;
    theta(l) = pi*rand(1)-pi/2;  %����·������
    
    % generate path directions
    flag = 1;
    while (flag)
        for pp = 1:l-1
            if (sqrt((psi(l)-psi(pp))^2*sys.N1*sys.N1+(theta(l)-theta(pp))^2*sys.N2*sys.N2)<1)  %���ƽǶȼ��
                psi(l) = pi*rand(1)-pi/2;
                theta(l) = pi*rand(1)-pi/2;  %��������·������
                break;
            end    
        end
        flag = 0;
    end
    
    % generate the path gains
    if (l==1)
        alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2)*sqrt(sys.kappa);  %��1��LoS·���Ĺ�һ������
    else
        alpha(l) = (normrnd(0, 1) + 1i*normrnd(0, 1))/sqrt(2);  %��2 -- L��NLoS·���Ĺ�һ������
    end
    operator = 1/sqrt(L-1+sys.kappa);  %����·����������
    
    
    h = h + 1/sqrt(sys.N)*operator*alpha(l)*kron(exp(1i*pi*sin(psi(l))*(0:1:(sys.N1-1)))', exp(1i*pi*cos(psi(l))*sin(theta(l))*(0:1:(sys.N2-1)))');
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

