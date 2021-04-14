clc;clear;close all;

%% ��ʼ��
system = SystemSet();  %ϵͳ������ʼ��

% ��ʼ��������Ͼ���
W = (normrnd(0,1,system.N,system.M)+1i*normrnd(0,1,system.N,system.M))/sqrt(2);  %������Ͼ��� N*M

% ��ʼ���ŵ�
[h, psi, theta, alpha] = RandomChannel_UPA(system);  %�ŵ���ģ
h_norm = norm(h);  %�ŵ�������ŷʽ����

% ��ʼ����λƫ��
phase = rand(1,system.N);  %����0-2pi�Ͼ��ȷֲ��������λƫ�� B*1
% phase = zeros(system.M,1);  %��ʼ����λƫ�� M*1
% for block_idx = 1:system.N_block
%     phase((block_idx-1)*system.N_RF+1:block_idx*system.N_RF,1) = phase_l(block_idx)*ones(system.N_RF,1);
% end   %����λƫ�ƾ���ֿ鸳ֵ 1��N_RF��M
phase_shift = diag(phase);  %��λƫ�ƾ���

% ��ʼ������
noise_var = 10^(-system.SNR/10)/system.N;  %��������
n_antenna = (normrnd(0,1,system.N,system.M) + 1i*normrnd(0,1,system.N,system.M))*sqrt(noise_var/2);  %���������˹������
n = diag(W'*n_antenna);  %�������������� M*1

% �����ź�
C = W'*exp(2i*pi*phase_shift);

y_true = C*h;  %������յ�Ƶ M*1
y = C*h + n;  %������յ�Ƶ M*1

h0 = initialization(y,W,system,10);
 
[phase_es,h_es,succ,lambda] = IR(y,h0,phase,phase,W,system);