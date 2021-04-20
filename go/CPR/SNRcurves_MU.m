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

SNR_list = -30:5:0;
SNR_list_size = numel(SNR_list);
M_default = 256;
N_RF_default = 16;
N_TS = 16;
sys_default = SystemSettings(M_default, N_RF_default, SNR_list(1));
options_default = MyDefaultOptions(sys_default);

K = N_RF_default;

%% the containers to put the results
h_norm = zeros(sys_default.N_trial, K);

SE_CRAF = zeros(SNR_list_size, sys_default.N_trial, K);
SNRLoss_CRAF = zeros(SNR_list_size, sys_default.N_trial, K);
success_rate_CRAF = zeros(SNR_list_size, 1);
SumRate_CRAF = zeros(SNR_list_size, sys_default.N_trial);

SE_PCCPR = zeros(SNR_list_size, sys_default.N_trial, sys_default.T+1, K);
SNRLoss_PCCPR = zeros(SNR_list_size, sys_default.N_trial, K);
success_rate_PCCPR = zeros(SNR_list_size, 1);
SumRate_PCCPR = zeros(SNR_list_size, sys_default.N_trial);

SE_PCCPR_SR = zeros(SNR_list_size, sys_default.N_trial, sys_default.T2+1, K);
SNRLoss_PCCPR_SR = zeros(SNR_list_size, sys_default.N_trial, K);
success_rate_PCCPR_SR = zeros(SNR_list_size, 1);
SumRate_PCCPR_SR = zeros(SNR_list_size, sys_default.N_trial);

SE_NOMP = zeros(SNR_list_size, sys_default.N_trial, K);
SNRLoss_NOMP = zeros(SNR_list_size, sys_default.N_trial, K);
success_rate_NOMP = zeros(SNR_list_size, 1);
SumRate_NOMP = zeros(SNR_list_size, sys_default.N_trial);

%% Start simulations
%% try Monte-Carlo simulations for N_trial instances
for trial_idx = 1:sys_default.N_trial
    %% set up a compressive phase retrieval problem
    % generate the precoding matrix
    if (strcmp(sys_default.array_type, 'UPA'))
        W = exp(1i*2*pi*rand(sys_default.N,sys_default.M));
        W_Gaussian = (normrnd(0,1,sys_default.N,sys_default.M)+1i*normrnd(0,1,sys_default.N,sys_default.M))/sqrt(2);
    elseif (strcmp(sys_default.array_type, 'ULA'))
        R = 4;
        Nb = M*R*R/N;
        [W, P] = (MultiArmMeasurements(N, R, Nb));
    end
    
    % generate the channel
    H = zeros(sys_default.N,K);
    for k_idx = 1:K
%       [h, theta1, theta2, alpha] = RandomChannel_UPA(sys_default);
        [h, theta1, theta2, alpha] = Channel_ULA(sys_default.N1, sys_default.N2, sys_default.L);
        H(:,k_idx) = h;
        h_norm(trial_idx, k_idx) = norm(h);
    end
    
    % generate random structured phase offsets
    phase_l = 2*pi*rand(options_default.N_block, K);
    phase = zeros(sys_default.M, K);
    for k_idx = 1:K
        for block_idx = 1:options_default.N_block
            phase((block_idx-1)*sys_default.block_size+1:block_idx*sys_default.block_size,k_idx) = phase_l(block_idx, k_idx)*ones(sys_default.block_size, 1);
        end
    end
    phase_noncoherent = 2*pi*rand(sys_default.M, K);
    
    for SNR_idx = 1:SNR_list_size       % go through all SNR's
        SNR = SNR_list(SNR_idx);        % inside the iteration, the number of measurement M is fixed
        
        % display the simulation progress
        disp(['trial=', num2str(trial_idx), '   SNR=', num2str(SNR)]);
        
        %% system and algorithm parameters
        sys = SystemSettings(M_default, N_RF_default, SNR);
        options = MyDefaultOptions(sys);% default options, refer to the file for details
        
        % generate noise
        noise_var_signal = 10^(-SNR/10);
        noise_var_pilot = 10^(-SNR/10)/N_TS;
        n = (normrnd(0, 1, sys_default.M, K) + 1i*normrnd(0, 1, sys_default.M, K))*sqrt(noise_var_pilot/2);

%         %% conventional CRAF with no joint phase
%         % calculate the received signal
%         H_CRAF = zeros(sys.N, K);
%         for k_idx = 1:K
%             h = H(:,k_idx);
%             y = (W_Gaussian'*h).*exp(1i*phase_noncoherent(:,k_idx)) + n(:,k_idx);
%             
%             [h_CRAF] = CRAF_CE(y, W_Gaussian, sys, options, h);
%             H_CRAF(:,k_idx) = h_CRAF;
%             SE = SE_rotate(h_CRAF, h);
%             SE_CRAF(SNR_idx, trial_idx, k_idx) = SE;
%             SNRLoss_CRAF(SNR_idx, trial_idx, k_idx) = abs((sign(h_CRAF)'*h)/(sign(h)'*h));
%         end
%         SumRate = EvaluateSumRate(sys,H,H_CRAF, noise_var_signal);
%         SumRate_CRAF(SNR_idx, trial_idx) = SumRate;
        
        %% partially coherent compressive channel estimation
        % calculate the received signal
%         H_PCCPR = zeros(sys.N, K);
%         for k_idx = 1:K
%             h = H(:,k_idx);
%             y = (W'*h).*exp(1i*phase(:,k_idx)) + n(:,k_idx);
% 
%             h0 = InitGuess_UPA(y, W, sys, options);
%             [h_PCCPR_list] = PCCPR_UPA(y, W, h0, sys, options);
%             H_PCCPR(:,k_idx) = h_PCCPR_list(:,sys_default.T+1);
%             SE_list = SE_rotate(h_PCCPR_list, h);
%             SE_PCCPR(SNR_idx, trial_idx, :, k_idx) = SE_list;
%             SNRLoss_PCCPR(SNR_idx, trial_idx, k_idx) = abs((sign(h_PCCPR_list(:,sys.T))'*h)/(sign(h)'*h));
%         end
%         SumRate = EvaluateSumRate(sys,H,H_PCCPR, noise_var_signal);
%         SumRate_PCCPR(SNR_idx, trial_idx) = SumRate;
        
        H_PCCPR_SR = zeros(sys.N, K);
        for k_idx = 1:K
            h = H(:,k_idx);
            y = (W'*h).*exp(1i*phase(:,k_idx)) + n(:,k_idx);

            h0 = InitGuess_UPA(y, W, sys, options);
            [h_PCCPR_SR_list] = PCCPR_SR_UPA(y, W, h0, sys, options);
            H_PCCPR_SR(:,k_idx) = h_PCCPR_SR_list(:,sys_default.T2+1);
            SE_list = SE_rotate(h_PCCPR_SR_list, h);
            SE_PCCPR_SR(SNR_idx, trial_idx, :, k_idx) = SE_list;
            SNRLoss_PCCPR_SR(SNR_idx, trial_idx, k_idx) = abs((sign(h_PCCPR_SR_list(:,sys.T))'*h)/(sign(h)'*h));
        end
        SumRate = EvaluateSumRate(sys,H,H_PCCPR_SR, noise_var_signal);
        SumRate_PCCPR_SR(SNR_idx, trial_idx) = SumRate;
        
        
        %% NOMP with perfect phase information
%         H_NOMP = zeros(sys.N, K);
%         for k_idx = 1:K
%             h = H(:,k_idx);
%             y = W'*h + n(:,k_idx);
%             
%             h_est_NOMP = NOMP_2D('2nd',sys_default.N1, sys_default.N2, W, y, options_default.Np, 1e-5, 5);
%             H_NOMP(:,k_idx) = h_est_NOMP;
%             SE_NOMP(SNR_idx, trial_idx, k_idx) = norm(h_est_NOMP - h, 2)^2;
%             SNRLoss_NOMP(SNR_idx, trial_idx, k_idx) = abs((sign(h_est_NOMP)'*h)/(sign(h)'*h));
%         end
%         SumRate = EvaluateSumRate(sys,H,H_NOMP, noise_var_signal);
%         SumRate_NOMP(SNR_idx, trial_idx) = SumRate;
    end
end

%% calculate the success rates
for SNR_idx = 1:SNR_list_size
%     success_rate_CRAF(SNR_idx) = mean(mean(SNRLoss_CRAF(SNR_idx,:,:)>0.5));
%     success_rate_PCCPR(SNR_idx) = mean(mean(SNRLoss_PCCPR(SNR_idx,:,:)>0.5));
    success_rate_PCCPR_SR(SNR_idx) = mean(mean(SNRLoss_PCCPR_SR(SNR_idx,:,:)>0.5));
%     success_rate_NOMP(SNR_idx) = mean(mean(SNRLoss_NOMP(SNR_idx,:,:)>0.5));
end

%% calculate the NMSE
% NMSE_CRAF = sum(sum(SE_CRAF,2),3)/sum(sum(h_norm.^2));
% NMSE_PCCPR = sum(sum(squeeze(SE_PCCPR(:,:,end,:)),2),3)/sum(sum(h_norm.^2));
NMSE_PCCPR_SR = sum(sum(squeeze(SE_PCCPR_SR(:,:,end,:)),2),3)/sum(sum(h_norm.^2));
% NMSE_NOMP = sum(sum(SE_NOMP,2),3)/sum(sum(h_norm.^2));


% figure;
% % plot(SNR_list,success_rate_PCCPR,'b-o','LineWidth',1.5);
% % hold on;
% plot(SNR_list,success_rate_PCCPR_SR,'c-*','LineWidth',1.5);
% hold on;
% plot(SNR_list,success_rate_NOMP(:),'g-v','LineWidth',1.5);
% hold on;
% plot(SNR_list,success_rate_CRAF(:),'r-^','LineWidth',1.5);
% legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
% grid on;
% xlabel('SNR');
% ylabel('Channel estimation success rate');

% figure;
% plot(SNR_list,mean(mean(SNRLoss_PCCPR, 2),3),'b-o','LineWidth',1.5);
% hold on;
% plot(SNR_list,mean(mean(SNRLoss_PCCPR_SR, 2),3),'c-*','LineWidth',1.5);
% hold on;
% plot(SNR_list,mean(mean(SNRLoss_NOMP, 2),3),'g-v','LineWidth',1.5);
% hold on;
% plot(SNR_list,mean(mean(SNRLoss_CRAF, 2),3),'r-^','LineWidth',1.5);
% legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
% grid on;
% xlabel('SNR');
% ylabel('mean SNR loss');

figure;
% semilogy(SNR_list,NMSE_PCCPR,'b-o','LineWidth',1.5);
% hold on;
semilogy(SNR_list,NMSE_PCCPR_SR,'c-*','LineWidth',1.5);
hold on;
% semilogy(SNR_list,NMSE_NOMP,'g-v','LineWidth',1.5);
% hold on;
% semilogy(SNR_list,NMSE_CRAF,'r-^','LineWidth',1.5);
% legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
grid on;
xlabel('SNR');
ylabel('NMSE performance');

% figure;
% plot(SNR_list,mean(SumRate_PCCPR,2),'b-o','LineWidth',1.5);
% hold on;
% plot(SNR_list,mean(SumRate_PCCPR_SR,2),'c-*','LineWidth',1.5);
% hold on;
% plot(SNR_list,mean(SumRate_NOMP,2),'g-v','LineWidth',1.5);
% hold on;
% plot(SNR_list,mean(SumRate_CRAF,2),'r-^','LineWidth',1.5);
% 
% legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
% grid on;
% xlabel('SNR');
% ylabel('NMSE performance');