clc;clear;close all;

%% 初始化
system = SystemSet();  %系统参数初始化

% 初始化编码组合矩阵
W = (normrnd(0,1,system.N,system.M)+1i*normrnd(0,1,system.N,system.M))/sqrt(2);  %编码组合矩阵 N*M

% 初始化信道
[h, psi, theta, alpha] = RandomChannel_UPA(system);  %信道建模
h_norm = norm(h);  %信道向量的欧式范数

% 初始化相位偏移
phase = rand(1,system.N);  %生成0-2pi上均匀分布的随机相位偏移 B*1
% phase = zeros(system.M,1);  %初始化相位偏移 M*1
% for block_idx = 1:system.N_block
%     phase((block_idx-1)*system.N_RF+1:block_idx*system.N_RF,1) = phase_l(block_idx)*ones(system.N_RF,1);
% end   %给相位偏移矩阵分块赋值 1：N_RF：M
phase_shift = diag(phase);  %相位偏移矩阵

% 初始化噪声
noise_var = 10^(-system.SNR/10)/system.N;  %噪声方差
n_antenna = (normrnd(0,1,system.N,system.M) + 1i*normrnd(0,1,system.N,system.M))*sqrt(noise_var/2);  %生成随机高斯复噪声
n = diag(W'*n_antenna);  %编码后的噪声向量 M*1

% 接收信号
C = W'*exp(2i*pi*phase_shift);

y_true = C*h;  %无噪接收导频 M*1
y = C*h + n;  %有噪接收导频 M*1

h0 = initialization(y,W,system,10);
 
[phase_es,h_es,succ,lambda] = IR(y,h0,phase,phase,W,system);