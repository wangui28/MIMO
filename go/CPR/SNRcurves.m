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

SNR_list = -5:2.5:15;
SNR_list_size = numel(SNR_list);
M_default = 256;
N_RF_default_4 = 4;
N_RF_default_16 = 16;
sys_default_RF4 = SystemSettings(M_default, N_RF_default_4, SNR_list(1));
options_default_RF4 = MyDefaultOptions(sys_default_RF4);
sys_default_RF16 = SystemSettings(M_default, N_RF_default_16, SNR_list(1));
options_default_RF16 = MyDefaultOptions(sys_default_RF16);


%% the containers to put the results
h_norm = zeros(sys_default_RF4.N_trial, 1);

SE_CRAF = zeros(SNR_list_size, sys_default_RF4.N_trial);
SNRLoss_CRAF = zeros(SNR_list_size, sys_default_RF4.N_trial);
success_rate_CRAF = zeros(SNR_list_size, 1);

SE_PCCPR_RF4 = zeros(SNR_list_size, sys_default_RF4.N_trial, sys_default_RF4.T+1);
SNRLoss_PCCPR_RF4 = zeros(SNR_list_size, sys_default_RF4.N_trial);
success_rate_PCCPR_RF4 = zeros(SNR_list_size, 1);

SE_PCCPR_SR_RF4 = zeros(SNR_list_size, sys_default_RF4.N_trial, sys_default_RF4.T2+1);
SNRLoss_PCCPR_SR_RF4 = zeros(SNR_list_size, sys_default_RF4.N_trial);
success_rate_PCCPR_SR_RF4 = zeros(SNR_list_size, 1);

SE_PCCPR_RF16 = zeros(SNR_list_size, sys_default_RF16.N_trial, sys_default_RF16.T+1);
SNRLoss_PCCPR_RF16 = zeros(SNR_list_size, sys_default_RF16.N_trial);
success_rate_PCCPR_RF16 = zeros(SNR_list_size, 1);

SE_PCCPR_SR_RF16 = zeros(SNR_list_size, sys_default_RF16.N_trial, sys_default_RF16.T2+1);
SNRLoss_PCCPR_SR_RF16 = zeros(SNR_list_size, sys_default_RF16.N_trial);
success_rate_PCCPR_SR_RF16 = zeros(SNR_list_size, 1);

SE_NOMP = zeros(SNR_list_size, sys_default_RF4.N_trial);
SNRLoss_NOMP = zeros(SNR_list_size, sys_default_RF4.N_trial);
success_rate_NOMP = zeros(SNR_list_size, 1);



%% Start simulations
%% try Monte-Carlo simulations for N_trial instances
for trial_idx = 1:sys_default_RF4.N_trial
    %% set up a compressive phase retrieval problem
    % generate the combining matrix
    if (strcmp(sys_default_RF4.array_type, 'UPA'))
        W = exp(1i*2*pi*rand(sys_default_RF4.N,sys_default_RF4.M));
        W_Gaussian = (normrnd(0,1,sys_default_RF4.N,sys_default_RF4.M)+1i*normrnd(0,1,sys_default_RF4.N,sys_default_RF4.M))/sqrt(2);
    elseif (strcmp(sys_default_RF4.array_type, 'ULA'))
        R = 4;
        Nb = M*R*R/N;
        [W, P] = (MultiArmMeasurements(N, R, Nb));
    end
        
    % generate the channel together with its parameters
    [h, theta1, theta2, alpha] = RandomChannel_UPA(sys_default_RF4);
    h_norm(trial_idx) = norm(h);
    
    % generate random structured phase offsets
    phase_l_RF4 = 2*pi*rand(options_default_RF4.N_block, 1);
    phase_RF4 = zeros(sys_default_RF4.M, 1);
    for block_idx = 1:options_default_RF4.N_block
        phase_RF4((block_idx-1)*sys_default_RF4.block_size+1:block_idx*sys_default_RF4.block_size,1) = phase_l_RF4(block_idx)*ones(sys_default_RF4.block_size, 1);
    end
    
    phase_l_RF16 = 2*pi*rand(options_default_RF16.N_block, 1);
    phase_RF16 = zeros(sys_default_RF16.M, 1);
    for block_idx = 1:options_default_RF16.N_block
        phase_RF16((block_idx-1)*sys_default_RF16.block_size+1:block_idx*sys_default_RF16.block_size,1) = phase_l_RF16(block_idx)*ones(sys_default_RF16.block_size, 1);
    end
    
    phase_noncoherent = 2*pi*rand(sys_default_RF4.M, 1);
    
    for SNR_idx = 1:SNR_list_size       % go through all SNR's
        SNR = SNR_list(SNR_idx);        % inside the iteration, the number of measurement M is fixed
        
        % display the simulation progress
        disp(['trial=', num2str(trial_idx), '   SNR=', num2str(SNR)]);
        
        %% system and algorithm parameters
        sys_1 = SystemSettings(M_default, N_RF_default_4, SNR);
        sys_2 = SystemSettings(M_default, N_RF_default_16, SNR);
        options_1 = MyDefaultOptions(sys_1);% default options, refer to the file for details
        options_2 = MyDefaultOptions(sys_2);
        
        % generate noise
        noise_var = 10^(-SNR/10);
        n = (normrnd(0, 1, sys_default_RF4.M, 1) + 1i*normrnd(0, 1, sys_default_RF4.M, 1))*sqrt(noise_var/2);

        %% conventional CRAF with no joint phase
        % calculate the received signal
        y = (W_Gaussian'*h).*exp(1i*phase_noncoherent) + n;
        
        [h_CRAF] = CRAF_CE(y, W_Gaussian, sys_1, options_1, h);
        SE = SE_rotate(h_CRAF, h);
        SE_CRAF(SNR_idx, trial_idx) = SE;
        SNRLoss_CRAF(SNR_idx, trial_idx) = abs((sign(h_CRAF)'*h)/(sign(h)'*h));
        
        %% partially coherent compressive channel estimation
            %% N_RF=4
        % calculate the received signal
        y = (W'*h).*exp(1i*phase_RF4) + n;
            
        h0 = InitGuess_UPA(y, W, sys_1, options_1);
        [h_PCCPR_list] = PCCPR_UPA(y, W, h0, sys_1, options_1);
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_RF4(SNR_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_RF4(SNR_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_1.T))'*h)/(sign(h)'*h));
        
        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y, W, h0, sys_1, options_1);
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_RF4(SNR_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_RF4(SNR_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_1.T2+1))'*h)/(sign(h)'*h));
        
            %% N_RF=16
        % calculate the received signal
        y = (W'*h).*exp(1i*phase_RF16) + n;
        
        % partially coherent compressive channel estimation
        h0 = InitGuess_UPA(y, W, sys_2, options_2);
        [h_PCCPR_list] = PCCPR_UPA(y, W, h0, sys_2, options_2);
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_RF16(SNR_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_RF16(SNR_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_2.T+1))'*h)/(sign(h)'*h));
        
        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y, W, h0, sys_2, options_2);
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_RF16(SNR_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_RF16(SNR_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_2.T2+1))'*h)/(sign(h)'*h));

        %% NOMP with perfect phase information
        y = W'*h + n;
        h_est_NOMP = NOMP_2D('2nd',sys_default_RF4.N1, sys_default_RF4.N2, W, y, options_default_RF4.Np, 1e-5, 5);
        SE_NOMP(SNR_idx, trial_idx) = norm(h_est_NOMP - h, 2)^2;
        SNRLoss_NOMP(SNR_idx, trial_idx) = abs((sign(h_est_NOMP)'*h)/(sign(h)'*h));
    end
    
end

%% calculate the success rates
for SNR_idx = 1:SNR_list_size
    success_rate_CRAF(SNR_idx) = mean(SNRLoss_CRAF(SNR_idx,:)>0.5);
    success_rate_PCCPR_RF4(SNR_idx) = mean(SNRLoss_PCCPR_RF4(SNR_idx,:)>0.5);
    success_rate_PCCPR_SR_RF4(SNR_idx) = mean(SNRLoss_PCCPR_SR_RF4(SNR_idx,:)>0.5);
    success_rate_PCCPR_RF16(SNR_idx) = mean(SNRLoss_PCCPR_RF16(SNR_idx,:)>0.5);
    success_rate_PCCPR_SR_RF16(SNR_idx) = mean(SNRLoss_PCCPR_SR_RF16(SNR_idx,:)>0.5);
    success_rate_NOMP(SNR_idx) = mean(SNRLoss_NOMP(SNR_idx,:)>0.5);
end

%% calculate the NMSE
NMSE_CRAF = sum(SE_CRAF,2)/sum(h_norm.^2);
NMSE_PCCPR_1 = sum(squeeze(SE_PCCPR_RF4(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_SR_1 = sum(squeeze(SE_PCCPR_SR_RF4(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_2 = sum(squeeze(SE_PCCPR_RF16(:,:,end)),2)/sum(h_norm.^2);
NMSE_PCCPR_SR_2 = sum(squeeze(SE_PCCPR_SR_RF16(:,:,end)),2)/sum(h_norm.^2);
NMSE_NOMP = sum(SE_NOMP,2)/sum(h_norm.^2);

save results.mat


figure;
plot(SNR_list,success_rate_PCCPR_RF4,'b-o','LineWidth',1.5);
hold on;
plot(SNR_list,success_rate_PCCPR_SR_RF4,'c-*','LineWidth',1.5);
hold on;
plot(SNR_list,success_rate_NOMP(:),'g-v','LineWidth',1.5);
hold on;
plot(SNR_list,success_rate_CRAF(:),'r-^','LineWidth',1.5);
hold on;
plot(SNR_list,success_rate_PCCPR_RF16,'b-o','LineWidth',1.5);
hold on;
plot(SNR_list,success_rate_PCCPR_SR_RF16,'c-*','LineWidth',1.5);
legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
grid on;
xlabel('SNR');
ylabel('Channel estimation success rate');

figure;
plot(SNR_list,mean(SNRLoss_PCCPR_RF4, 2),'b-o','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SNRLoss_PCCPR_SR_RF4, 2),'c-*','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SNRLoss_NOMP, 2),'g-v','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SNRLoss_CRAF, 2),'r-^','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SNRLoss_PCCPR_RF16, 2),'b-o','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SNRLoss_PCCPR_SR_RF16, 2),'c-*','LineWidth',1.5);
legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
grid on;
xlabel('SNR');
ylabel('mean SNR loss');

figure;
semilogy(SNR_list,NMSE_PCCPR_1,'b-o','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_PCCPR_SR_1,'c-*','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_NOMP,'g-v','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_CRAF,'r-^','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_PCCPR_2,'b-o','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_PCCPR_SR_2,'c-*','LineWidth',1.5);
legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF')
grid on;
xlabel('SNR');
ylabel('NMSE performance');