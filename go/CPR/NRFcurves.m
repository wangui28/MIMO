clear all;clc;close all;

%% overall settings, set up the system parameters
%       OPEN system_settings.m, CHECK OR MODIFY the system parameters 
%       before simulation


% struct "sys" includes the following parameters:
%
%    array parameters:
%       sys_default.array_type:     the antenna array type, typically "UPA"
%       sys_default.N1:             (array_type=="UPA") UPA size in azimuth
%       sys_default.N2:             (array_type=="UPA") UPA size in elevation
%       sys_default.N:              the number of antennas
%       sys_default.block_size      the size of the phase structure block, which 
%                           equals to the number of RF chains
%
%   channel parameters:
%       sys_default.L:              the number of paths (including LoS & NLoS)
%
%   simulation settings
%       sys_default.SNR:            SNR for test
%       sys_default.M_list:         the list of the numbers of measurements for
%                           which this simulation is going to test
%       sys_default.M_list_size:    the size of the above list
%       sys_default.N_trial;        number of Monte-Carlo trials
%

N_RF_list = [1,2,4,8,16];  %RF链列表
N_RF_list_size = numel(N_RF_list);  %列表中元素个数
M_default = 256;  %测量数M
SNR_default = 10;  %信噪比
sys_default = SystemSettings(M_default, N_RF_list(1), SNR_default);  %系统模型参数
options_default = MyDefaultOptions(sys_default);  %预置选项


%% the containers to put the results
h_norm = zeros(N_RF_list_size, sys_default.N_trial);
SE_CRAF = zeros(N_RF_list_size, sys_default.N_trial);
SNRLoss_CRAF = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_CRAF = zeros(N_RF_list_size, 1);
SE_PCCPR_1 = zeros(N_RF_list_size, sys_default.N_trial, sys_default.T+1);
SNRLoss_PCCPR_1 = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_PCCPR_1 = zeros(N_RF_list_size, 1);
SE_PCCPR_SR_1 = zeros(N_RF_list_size, sys_default.N_trial, sys_default.T2+1);
SNRLoss_PCCPR_SR_1 = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_PCCPR_SR_1 = zeros(N_RF_list_size, 1);
SE_PCCPR_2 = zeros(N_RF_list_size, sys_default.N_trial, sys_default.T+1);
SNRLoss_PCCPR_2 = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_PCCPR_2 = zeros(N_RF_list_size, 1);
SE_PCCPR_SR_2 = zeros(N_RF_list_size, sys_default.N_trial, sys_default.T2+1);
SNRLoss_PCCPR_SR_2 = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_PCCPR_SR_2 = zeros(N_RF_list_size, 1);
SE_OMP = zeros(N_RF_list_size, sys_default.N_trial);
SNRLoss_OMP = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_OMP = zeros(N_RF_list_size, 1);
SE_NOMP = zeros(N_RF_list_size, sys_default.N_trial);
SNRLoss_NOMP = zeros(N_RF_list_size, sys_default.N_trial);
success_rate_NOMP = zeros(N_RF_list_size, 1);
% Angle_err = zeros(SNR_list_size, sys_default.N_trial);
% SNRLoss_FBA = zeros(SNR_list_size, sys_default.N_trial);
% success_rate_FBA = zeros(SNR_list_size, 1);


%% Start simulations

for NRF_idx = 1:N_RF_list_size       % go through all N_RF's
    N_RF = N_RF_list(NRF_idx);        % inside the iteration, the number of measurement M is fixed
                                    
    %% system and algorithm parameters
    sys = SystemSettings(M_default, N_RF, SNR_default);
    options = MyDefaultOptions(sys);% default options, refer to the file for details

    
    %% try Monte-Carlo simulations for N_trial instances
    for trial_idx = 1:sys_default.N_trial
        
        %% set up a compressive phase retrieval problem
        % display the simulation progress
        disp(['N_RF=', num2str(N_RF), '  trial=', num2str(trial_idx)]);
        
        % generate the combining matrix
        if (strcmp(sys.array_type, 'UPA'))
            W = exp(1i*2*pi*rand(sys.N,sys.M));  %随机编码组合矩阵
            W_Gaussian = (normrnd(0,1,sys_default.N,sys_default.M)+1i*normrnd(0,1,sys_default.N,sys_default.M))/sqrt(2);  %服从高斯分布的编码组合矩阵
        elseif (strcmp(sys.array_type, 'ULA'))
            R = 4;
            Nb = M*R*R/N;
            [W, P] = (MultiArmMeasurements(N, R, Nb));
        end
        
        D1 = exp(-1i*2*pi*(0:1:(sys.N1-1))'*(0:1:(options.n1_dic-1))/options.n1_dic)/sqrt(sys.N1);
        D2 = exp(-1i*2*pi*(0:1:(sys.N2-1))'*(0:1:(options.n2_dic-1))/options.n2_dic)/sqrt(sys.N2);
        D = kron(D1, D2);  %构造离散域网络
        A = W'*D;  %离散组合矩阵
        
        % generate the channel together with its parameters
        [h, theta1, theta2, alpha] = RandomChannel_UPA(sys_default);  %%信道建模
        h_norm(NRF_idx, trial_idx) = norm(h);
        
        % generate random structured phase offsets
        phase_l = 2*pi*rand(options.N_block, 1);  %生成0-2pi上均匀分布的随机相位偏移
        phase = zeros(sys.M, 1);  %生成相位偏移矩阵 256*1
        for block_idx = 1:options.N_block
            phase((block_idx-1)*sys.block_size+1:block_idx*sys.block_size,1) = phase_l(block_idx)*ones(sys.block_size, 1);
        end  %给相位偏移矩阵分块赋值 1：Nrf：256
        
        % generate noise
        noise_var = 10^(-SNR_default/10);
        n = (normrnd(0, 1, sys.M, 1) + 1i*normrnd(0, 1, sys.M, 1))*sqrt(noise_var/2);  %生成随机高斯复噪声
        
        % calculate the received signal
        y = (W'*h).*exp(1i*phase) + n;  %随机
        y_noncoherent = (W_Gaussian'*h).*exp(1i*2*pi*rand(sys_default.M, 1)) + n;  %高斯分布

        %% conventional CRAF with no joint phase
        if (NRF_idx==1)
            [h_CRAF] = CRAF_CE(y_noncoherent, W_Gaussian, sys, options, h);  %CRAF
            SE = SE_rotate(h_CRAF, h);
            SE_CRAF(NRF_idx, trial_idx) = SE;
            SNRLoss_CRAF(NRF_idx, trial_idx) = abs((sign(h_CRAF)'*h)/(sign(h)'*h));
        else
            SE_CRAF(NRF_idx, trial_idx) = SE_CRAF(1, trial_idx);
            SNRLoss_CRAF(NRF_idx, trial_idx) = SNRLoss_CRAF(1, trial_idx);
        end
        
        %% partially coherent compressive channel estimation
        h0 = InitGuess_UPA(y, W, sys, options);
        [h_PCCPR_list] = PCCPR_UPA(y, W, h0, sys, options);  %on-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_1(NRF_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_1(NRF_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys.T+1))'*h)/(sign(h)'*h));
        
        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y, W, h0, sys, options);  %off-grid PC-CPR
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_1(NRF_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_1(NRF_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys.T2+1))'*h)/(sign(h)'*h));

        %% OMP with perfect phase information
        if (NRF_idx == 1)
            y = W'*h + n;
            h_est_OMP = D*OMP(A, y, options.k, 1e-4);  %OMP
            SE_OMP(NRF_idx, trial_idx) = norm(h_est_OMP - h, 2)^2;
            SNRLoss_OMP(NRF_idx, trial_idx) = abs((sign(h_est_OMP)'*h)/(sign(h)'*h));
        else
            SE_OMP(NRF_idx, trial_idx) = SE_OMP(1, trial_idx);
            SNRLoss_OMP(NRF_idx, trial_idx) = SNRLoss_OMP(1, trial_idx);
        end
        
        %% NOMP with perfect phase information
        if (NRF_idx == 1)
            y = W'*h + n;
            h_est_NOMP = NOMP_2D('2nd',sys.N1, sys.N2, W, y, options_default.Np, 1e-5, 5);  %NOMP
            SE_NOMP(NRF_idx, trial_idx) = norm(h_est_NOMP - h, 2)^2;
            SNRLoss_NOMP(NRF_idx, trial_idx) = abs((sign(h_est_NOMP)'*h)/(sign(h)'*h));
        else
            SE_NOMP(NRF_idx, trial_idx) = SE_NOMP(1, trial_idx);
            SNRLoss_NOMP(NRF_idx, trial_idx) = SNRLoss_NOMP(1, trial_idx);
        end
        
          %% Fast MmWave Beam Alignment
%         TT = FBA(y, W, P, R, N, Nb);
%         [~, TTi] = max(TT);
%         Angle_err(M_iter, l) = (mod(TTi/N/2-0.5,1)-0.5)*2*pi - theta(1);
%         SNRLoss_FBA(M_iter, l) = abs((exp(1i*(mod(TTi/N/2-0.5,1)-0.5)*2*pi*(0:1:(N-1))'))'*h/(sign(h)'*h));
    end
    
    %% calculate the success rates
    success_rate_CRAF(NRF_idx) = mean(SNRLoss_CRAF(NRF_idx,:)>0.5);
    success_rate_PCCPR_1(NRF_idx) = mean(SNRLoss_PCCPR_1(NRF_idx,:)>0.5);
    success_rate_PCCPR_SR_1(NRF_idx) = mean(SNRLoss_PCCPR_SR_1(NRF_idx,:)>0.5);
    success_rate_NOMP(NRF_idx) = mean(SNRLoss_NOMP(NRF_idx,:)>0.5);
    success_rate_OMP(NRF_idx) = mean(SNRLoss_OMP(NRF_idx,:)>0.5);
end

%% calculate the NMSE
NMSE_CRAF = sum(SE_CRAF,2)/sum(h_norm(1,:).^2,2);
NMSE_PCCPR_1 = sum(squeeze(SE_PCCPR_1(:,:,end)),2)./sum(h_norm.^2,2);
NMSE_PCCPR_SR_1 = sum(squeeze(SE_PCCPR_SR_1(:,:,end)),2)./sum(h_norm.^2,2);
NMSE_OMP = sum(SE_OMP,2)/sum(h_norm(1,:).^2,2);
NMSE_NOMP = sum(SE_NOMP,2)/sum(h_norm(1,:).^2,2);

save results.mat


figure;
plot(N_RF_list,success_rate_PCCPR_1,'b-o','LineWidth',1.5);
hold on;
plot(N_RF_list,success_rate_PCCPR_SR_1,'c-*','LineWidth',1.5);
hold on;
plot(N_RF_list,success_rate_CRAF(:),'r-^','LineWidth',1.5);
hold on;
plot(N_RF_list,success_rate_OMP(:),'m-x','LineWidth',1.5,'Markersize',10);
hold on;
plot(N_RF_list,success_rate_NOMP(:),'g-v','LineWidth',1.5);

% hold on;
% plot(NRF_list,success_rate_PCCPR_2,'b-o','LineWidth',1.5);
% hold on;
% plot(NRF_list,success_rate_PCCPR_SR_2,'c-*','LineWidth',1.5);
legend('proposed PC-CPR (on-grid)','proposed PC-CPR (off-grid)','CRAF [17]','OMP assuming known phase','NOMP assuming known phase [30]')
grid on;
xlim([0 18]);
xlabel('{\it N}_{\rm RF}');
ylabel('Channel estimation success rate');

figure;
plot(N_RF_list,median(-10*log10(SNRLoss_PCCPR_1), 2),'b-o','LineWidth',1.5);
hold on;
plot(N_RF_list,median(-10*log10(SNRLoss_PCCPR_SR_1), 2),'c-*','LineWidth',1.5);
hold on;
plot(N_RF_list,median(-10*log10(SNRLoss_CRAF), 2),'r-^','LineWidth',1.5);
hold on;
plot(N_RF_list,median(-10*log10(SNRLoss_OMP), 2),'m-x','LineWidth',1.5,'Markersize',10);
hold on;
plot(N_RF_list,median(-10*log10(SNRLoss_NOMP), 2),'g-v','LineWidth',1.5);
% hold on;
% plot(NRF_list,mean(SNRLoss_PCCPR_2, 2),'b-o','LineWidth',1.5);
% hold on;
% plot(NRF_list,mean(SNRLoss_PCCPR_SR_2, 2),'c-*','LineWidth',1.5);
legend('proposed PC-CPR (on-grid)','proposed PC-CPR (off-grid)','CRAF [17]','OMP assuming known phase','NOMP assuming known phase [30]')
grid on;
xlim([0 18]);
xlabel('{\it N}_{\rm RF}');
ylabel('median SNR loss (dB)');

figure;
semilogy(N_RF_list,NMSE_PCCPR_1,'b-o','LineWidth',1.5);
hold on;
semilogy(N_RF_list,NMSE_PCCPR_SR_1,'c-*','LineWidth',1.5);
hold on;
semilogy(N_RF_list,NMSE_CRAF,'r-^','LineWidth',1.5);
hold on;
semilogy(N_RF_list,NMSE_OMP,'m-x','LineWidth',1.5,'Markersize',10);
hold on;
semilogy(N_RF_list,NMSE_NOMP,'g-v','LineWidth',1.5);
% hold on;
% semilogy(NRF_list,NMSE_PCCPR_2,'b-o','LineWidth',1.5);
% hold on;
% semilogy(NRF_list,NMSE_PCCPR_SR_2,'c-*','LineWidth',1.5);
legend('proposed PC-CPR (on-grid)','proposed PC-CPR (off-grid)','CRAF [17]','OMP assuming known phase','NOMP assuming known phase [30]')
grid on;
xlim([0 18]);
ylim([0.001 1]);
xlabel('{\it N}_{\rm RF}');
ylabel('NMSE performance');