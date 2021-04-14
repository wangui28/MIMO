clear all;clc;close all;

%% overall settings, set up the system parameters
%       OPEN system_settings.m, CHECK OR MODIFY the system parameters
%       before simulation


% struct  "sys" includes the following parameters:
%
%    array parameters:
%       sys.array_type:     the antenna array type, typically "UPA"
%       sys.N1:             (array_type=="UPA") UPA size in azimuth
%       sys.N2:             (array_type=="UPA") UPA size in elevation
%       sys.N:              the number of antennas
%       sys.block_size      the size of the phase structure block, which
%                           equals to the number of RF chains
%
%   channel parameters:
%       sys.L:              the number of paths (including LoS & NLoS)
%
%   simulation settings
%       sys.SNR:            SNR for test
%       sys.M_list:         the list of the numbers of measurements for
%                           which this simulation is going to test
%       sys.M_list_size:    the size of the above list
%       sys.N_trial;        number of Monte-Carlo trials
%

M_list = (32:32:256);  %测试数量列表
M_list_size = numel(M_list);  %列表中元素个数
SNR_default = 10;  %预置信噪比SNR
N_RF_default = 16;   %预置RF链数量
sys_default = SystemSettings(max(M_list), N_RF_default, SNR_default);  %系统模型参数
options_default = MyDefaultOptions(sys_default);  %预置选项
%% the containers of the results
h_norm = zeros(sys_default.N_trial,1);  %生成结果存储矩阵

%on-grid PC-CPR非相干预处理
SE_PCCPR_bm = zeros(M_list_size, sys_default.N_trial, sys_default.T+1);  %生成均方误差存储矩阵
SNRLoss_PCCPR_bm = zeros(M_list_size, sys_default.N_trial);  %生成损失信噪比存储矩阵
success_rate_PCCPR_bm = zeros(M_list_size, 1);  %生成成功率存储矩阵

%off-grid PC-CPR非相干预处理
SE_PCCPR_SR_bm = zeros(M_list_size, sys_default.N_trial, sys_default.T2+1);
SNRLoss_PCCPR_SR_bm = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_SR_bm = zeros(M_list_size, 1);

%on-grid PC-CPR最佳预处理
SE_PCCPR_opt = zeros(M_list_size, sys_default.N_trial, sys_default.T+1);
SNRLoss_PCCPR_opt = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_opt = zeros(M_list_size, 1);

%off-grid PC-CPR最佳预处理
SE_PCCPR_SR_opt = zeros(M_list_size, sys_default.N_trial, sys_default.T2+1);
SNRLoss_PCCPR_SR_opt = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_SR_opt = zeros(M_list_size, 1);

%on-grid PC-CPR MM预处理
SE_PCCPR_MM = zeros(M_list_size, sys_default.N_trial, sys_default.T+1);
SNRLoss_PCCPR_MM = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_MM = zeros(M_list_size, 1);

%off-grid PC-CPR  MM预处理
SE_PCCPR_SR_MM = zeros(M_list_size, sys_default.N_trial, sys_default.T2+1);
SNRLoss_PCCPR_SR_MM = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_SR_MM = zeros(M_list_size, 1);

%on-grid PC-CPR随机预处理
SE_PCCPR_rnd = zeros(M_list_size, sys_default.N_trial, sys_default.T+1);
SNRLoss_PCCPR_rnd = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_rnd = zeros(M_list_size, 1);

%off-grid PC-CPR随机预处理
SE_PCCPR_SR_rnd = zeros(M_list_size, sys_default.N_trial, sys_default.T2+1);
SNRLoss_PCCPR_SR_rnd = zeros(M_list_size, sys_default.N_trial);
success_rate_PCCPR_SR_rnd = zeros(M_list_size, 1);


%% Start simulations
%% try Monte-Carlo simulations for N_trial instances
for trial_idx = 1:sys_default.N_trial
    %% set up a compressive phase retrieval problem
    % generate the combining matrix
    if (strcmp(sys_default.array_type, 'UPA'))  %均匀平面阵列
% 	    W = exp(1i*2*pi*rand(sys_default.N,sys_default.M));
        W = (normrnd(0,1,sys_default.N,sys_default.M)+1i*normrnd(0,1,sys_default.N,sys_default.M))/sqrt(2);  %生成编码组合矩阵
    elseif (strcmp(sys_default.array_type, 'ULA'))  %均匀线性阵列
        R = 4;
        Nb = sys_default.M*R*R/sys_default.N; 
        [W, P] = (MultiArmMeasurements(sys_default.N, R, Nb));
    end

    % generate the channel together with its parameters
    [h, psi, theta, alpha] = RandomChannel_UPA(sys_default);  %信道建模
    h_norm(trial_idx) = norm(h);  %信道向量的欧式范数

    % generate random structured phase offsets
    phase_l = 2*pi*rand(options_default.N_block, 1);  %生成0-2pi上均匀分布的随机相位偏移 B*1
    phase = zeros(sys_default.M, 1);  %生成相位偏移矩阵 Max*1
    for block_idx = 1:options_default.N_block
        phase((block_idx-1)*sys_default.block_size+1:block_idx*sys_default.block_size,1) = phase_l(block_idx)*ones(sys_default.block_size, 1);
    end   %给相位偏移矩阵分块赋值 1：Nrf：Max

    % generate noise
    noise_var = 10^(-sys_default.SNR/10)/sys_default.N;  %噪声方差
    n_antenna = (normrnd(0, 1, sys_default.N, sys_default.M) + 1i*normrnd(0, 1, sys_default.N, sys_default.M))*sqrt(noise_var/2);  %生成随机高斯复噪声
    n = diag(W'*n_antenna);  %编码后的噪声向量 Max*1

    % calculate the received signal
    Y_nf = (W'*h).*exp(1i*phase);  %无噪接收导频 Max*1
    y = (W'*h).*exp(1i*phase) + n;  %有噪接收导频 Max*1

    for M_idx = 1:M_list_size       % go through all M's
        M = M_list(M_idx);          % inside the iteration, the number of measurement M is fixed
                                    %取不同固定M运行
        %% display the simulation progress
        disp(['M=', num2str(M), '  trial=', num2str(trial_idx)]);  %窗口显示仿真进度

        %% system and algorithm parameters
        sys_M = SystemSettings(M, N_RF_default, SNR_default);  %仿真模型参数
        options_M = MyDefaultOptions(sys_M);  %算法预置参数

        y_M = y(1:M);
        y_nf = Y_nf(1:M);
        W_M = W(:,1:M);
        n_M = n(1:M);  %取对应M个测量值

        %% partially coherent compressive channel estimation
        %% non-coherent preprocessing function非相干预处理
        sys_M_benchmark = SystemSettings(M, 1, SNR_default);  %Nrf=1时，系统模型为非相干
        options_M_benchmark = MyDefaultOptions(sys_M_benchmark);
        h0 = InitGuess_UPA(y_M, W_M, sys_M_benchmark, options_M_benchmark, sys_M.N*noise_var);  %初始化角度域信道向量
%         h0 = InitGuess_1(y_M*sqrt(M), W_M*sqrt(M), sys_M, options_M, noise_var/2*M);
%         h0 = InitGuess_1(y_M, W_M, sys_M, options_M, sys_M.N*noise_var, y_nf, n_M);
        [h_PCCPR_list] = PCCPR_UPA(y_M, W_M, h0, sys_M, options_M);  %on-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_list, h);  
        SE_PCCPR_bm(M_idx, trial_idx, :) = SE_list;  %计算均方误差
        SNRLoss_PCCPR_bm(M_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_M.T))'*h)/(sign(h)'*h));  %计算损失信噪比

        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y_M, W_M, h0, sys_M, options_M);  %off-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_bm(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_bm(M_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_M.T2))'*h)/(sign(h)'*h));
        
        %% partially coherent compressive channel estimation
        %% optimal preprocessing function最佳预处理
        h0 = InitGuess_UPA(y_M, W_M, sys_M, options_M, sys_M.N*noise_var);
%         h0 = InitGuess_1(y_M*sqrt(M), W_M*sqrt(M), sys_M, options_M, noise_var/2*M);
%         h0 = InitGuess_1(y_M, W_M, sys_M, options_M, sys_M.N*noise_var, y_nf, n_M);
        [h_PCCPR_list] = PCCPR_UPA(y_M, W_M, h0, sys_M, options_M);  %on-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_opt(M_idx, trial_idx, :) = SE_list;  %计算均方误差
        SNRLoss_PCCPR_opt(M_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_M.T))'*h)/(sign(h)'*h));  %计算损失信噪比

        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y_M, W_M, h0, sys_M, options_M);  %off-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_opt(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_opt(M_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_M.T2))'*h)/(sign(h)'*h));

        %% MM preprocessing function Mondelli&Montanari的预处理
        h0 = InitGuess_1(y_M, W_M, sys_M, options_M, sys_M.N*noise_var);
%         h0 = InitGuess_1(y_M*sqrt(M), W_M*sqrt(M), sys_M, options_M, noise_var/2*M);
%         h0 = InitGuess_1(y_M, W_M, sys_M, options_M, sys_M.N*noise_var, y_nf, n_M);
        [h_PCCPR_list] = PCCPR_UPA(y_M, W_M, h0, sys_M, options_M);  %on-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_MM(M_idx, trial_idx, :) = SE_list;  %计算均方误差
        SNRLoss_PCCPR_MM(M_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_M.T))'*h)/(sign(h)'*h));  %计算损失信噪比

        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y_M, W_M, h0, sys_M, options_M);  %off-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_MM(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_MM(M_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_M.T2))'*h)/(sign(h)'*h));
        
        %% random initialization随机初始化
        h0 = (normrnd(0,1,sys_M.N,1)+1i*normrnd(0,1,sys_M.N,1))/sqrt(sys_M.N*2);
        [h_PCCPR_list] = PCCPR_UPA(y_M, W_M, h0, sys_M, options_M);  %on-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_rnd(M_idx, trial_idx, :) = SE_list;  %计算均方误差
        SNRLoss_PCCPR_rnd(M_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_M.T))'*h)/(sign(h)'*h));  %计算损失信噪比

%         [h_PCCPR_SR_list] = PCCPR_SR_UPA(y_M, W_M, h0, sys_M, options_M);
%         SE_list = SE_rotate(h_PCCPR_SR_list, h);
%         SE_PCCPR_SR_rnd(M_idx, trial_idx, :) = SE_list;
%         SNRLoss_PCCPR_SR_rnd(M_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_M.T2))'*h)/(sign(h)'*h));
    end
end

for M_idx = 1:M_list_size
    %% calculate the success rates计算成功率
    success_rate_PCCPR_bm(M_idx) = mean(SNRLoss_PCCPR_bm(M_idx,:)>0.5);
    success_rate_PCCPR_SR_bm(M_idx) = mean(SNRLoss_PCCPR_SR_bm(M_idx,:)>0.5);
    success_rate_PCCPR_opt(M_idx) = mean(SNRLoss_PCCPR_opt(M_idx,:)>0.5);
    success_rate_PCCPR_SR_opt(M_idx) = mean(SNRLoss_PCCPR_SR_opt(M_idx,:)>0.5);
    success_rate_PCCPR_MM(M_idx) = mean(SNRLoss_PCCPR_MM(M_idx,:)>0.5);
    success_rate_PCCPR_SR_MM(M_idx) = mean(SNRLoss_PCCPR_SR_MM(M_idx,:)>0.5);
    success_rate_PCCPR_rnd(M_idx) = mean(SNRLoss_PCCPR_rnd(M_idx,:)>0.5);
    success_rate_PCCPR_SR_rnd(M_idx) = mean(SNRLoss_PCCPR_SR_rnd(M_idx,:)>0.5);
end

%% calculate the NMSE计算最小均方误差
NMSE_PCCPR_bm = sum(squeeze(SE_PCCPR_bm(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_SR_bm = sum(squeeze(SE_PCCPR_SR_bm(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_opt = sum(squeeze(SE_PCCPR_opt(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_SR_opt = sum(squeeze(SE_PCCPR_SR_opt(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_MM = sum(squeeze(SE_PCCPR_MM(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_SR_MM = sum(squeeze(SE_PCCPR_SR_MM(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_rnd = sum(squeeze(SE_PCCPR_rnd(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_SR_rnd = sum(squeeze(SE_PCCPR_SR_rnd(:,:,end)),2)/sum(h_norm.^2);
save results.mat


figure;
plot(M_list,success_rate_PCCPR_bm,'k-d','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_opt,'r-o','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_MM,'c-^','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_rnd,'b-v','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_bm,'k-x','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_opt,'r-x','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_MM,'c-x','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_rnd,'b-x','LineWidth',1.5);
legend('benchmark','opt','MM','rnd','benchmark','opt','MM','rnd');
grid on;
xlabel('M');
ylabel('Channel estimation success rate');

figure;
plot(M_list,mean(SNRLoss_PCCPR_bm, 2),'k-d','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_opt, 2),'r-o','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_MM, 2),'c-^','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_rnd, 2),'b-v','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_SR_bm, 2),'k-x','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_SR_opt, 2),'r-x','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_SR_MM, 2),'c-x','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_SR_rnd, 2),'b-x','LineWidth',1.5);
legend('benchmark','opt','MM','rnd','benchmark','opt','MM','rnd');
grid on;
xlabel('M');
ylabel('mean SNR loss');

figure;
semilogy(M_list,NMSE_PCCPR_bm,'k-d','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_opt,'r-o','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_MM,'c-^','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_rnd,'b-v','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_bm,'k-x','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_opt,'r-x','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_MM,'c-x','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_rnd,'b-x','LineWidth',1.5);
legend('benchmark','opt','MM','rnd','benchmark','opt','MM','rnd');
grid on;
xlabel('M');
ylabel('NMSE performance');