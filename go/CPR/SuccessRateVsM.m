clear all;clc;close all;

%% overall settings, set up the system parameters
%       OPEN system_settings.m, CHECK OR MODIFY the system parameters
%       before simulation


% struct "sys" includes the following parameters:
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

M_list = (32:32:256);
M_list_size = numel(M_list);
SNR_default = 10;
N_RF_default = 4;
sys_default = SystemSettings(max(M_list), N_RF_default, SNR_default);
options_default = MyDefaultOptions(sys_default);

N_RF_4 = 4;
sys_RF4 = SystemSettings(max(M_list), N_RF_4, SNR_default);
options_RF4 = MyDefaultOptions(sys_RF4);
N_RF_16 = 16;
sys_RF16 = SystemSettings(max(M_list), N_RF_16, SNR_default);
options_RF16 = MyDefaultOptions(sys_RF16);

%% the containers of the results
h_norm = zeros(sys_default.N_trial,1);

SE_CRAF = zeros(M_list_size, sys_default.N_trial);
SNRLoss_CRAF = zeros(M_list_size, sys_default.N_trial);
success_rate_CRAF = zeros(M_list_size, 1);

SE_PCCPR_RF4 = zeros(M_list_size, sys_RF4.N_trial, sys_RF4.T+1);
SNRLoss_PCCPR_RF4 = zeros(M_list_size, sys_RF4.N_trial);
success_rate_PCCPR_RF4 = zeros(M_list_size, 1);

SE_PCCPR_SR_RF4 = zeros(M_list_size, sys_RF4.N_trial, sys_RF4.T2+1);
SNRLoss_PCCPR_SR_RF4 = zeros(M_list_size, sys_RF4.N_trial);
success_rate_PCCPR_SR_RF4 = zeros(M_list_size, 1);

SE_PCCPR_RF16 = zeros(M_list_size, sys_RF16.N_trial, sys_RF16.T+1);
SNRLoss_PCCPR_RF16 = zeros(M_list_size, sys_RF16.N_trial);
success_rate_PCCPR_RF16 = zeros(M_list_size, 1);

SE_PCCPR_SR_RF16 = zeros(M_list_size, sys_RF16.N_trial, sys_RF16.T2+1);
SNRLoss_PCCPR_SR_RF16 = zeros(M_list_size, sys_RF16.N_trial);
success_rate_PCCPR_SR_RF16 = zeros(M_list_size, 1);

SE_NOMP = zeros(M_list_size, sys_default.N_trial);
SNRLoss_NOMP = zeros(M_list_size, sys_default.N_trial);
success_rate_NOMP = zeros(M_list_size, 1);


%% Start simulations
%% try Monte-Carlo simulations for N_trial instances
for trial_idx = 1:sys_default.N_trial
    %% set up a compressive phase retrieval problem
    % generate the combining matrix
    if (strcmp(sys_default.array_type, 'UPA'))
 	    W = exp(1i*2*pi*rand(sys_default.N,sys_default.M));
        W_Gaussian = (normrnd(0,1,sys_default.N,sys_default.M)+1i*normrnd(0,1,sys_default.N,sys_default.M))/sqrt(2);
    elseif (strcmp(sys_default.array_type, 'ULA'))
        R = 4;
        Nb = sys_default.M*R*R/sys_default.N;
        [W, P] = (MultiArmMeasurements(sys_default.N, R, Nb));
    end

    % generate the channel together with its parameters
    [h, theta1, theta2, alpha] = RandomChannel_UPA(sys_default);
    h_norm(trial_idx) = norm(h);

    % generate random structured phase offsets
    phase_l_RF4 = 2*pi*rand(options_RF4.N_block, 1);
    phase_RF4 = zeros(sys_RF4.M, 1);
    for block_idx = 1:options_RF4.N_block
        phase_RF4((block_idx-1)*sys_RF4.block_size+1:block_idx*sys_RF4.block_size,1) = phase_l_RF4(block_idx)*ones(sys_RF4.block_size, 1);
    end
    
    phase_l_RF16 = 2*pi*rand(options_RF16.N_block, 1);
    phase_RF16 = zeros(sys_RF16.M, 1);
    for block_idx = 1:options_RF16.N_block
        phase_RF16((block_idx-1)*sys_RF16.block_size+1:block_idx*sys_RF16.block_size,1) = phase_l_RF4(block_idx)*ones(sys_RF16.block_size, 1);
    end

    % generate noise
    noise_var = 10^(-sys_default.SNR/10);
    n = (normrnd(0, 1, sys_default.M, 1) + 1i*normrnd(0, 1, sys_default.M, 1))*sqrt(noise_var/2);

    % calculate the received signal
    Y_nf_RF4 = (W'*h).*exp(1i*phase_RF4);
    y_RF4 = (W'*h).*exp(1i*phase_RF4) + n;
    Y_nf_RF16 = (W'*h).*exp(1i*phase_RF16);
    y_RF16 = (W'*h).*exp(1i*phase_RF16) + n;
    y_noncoherent = (W_Gaussian'*h).*exp(1i*2*pi*rand(sys_default.M, 1)) + n;
    
    for M_idx = 1:M_list_size       % go through all M's
        M = M_list(M_idx);          % inside the iteration, the number of measurement M is fixed
        sys_M_RF4 = SystemSettings(M, N_RF_4, SNR_default);
        options_M_RF4 = MyDefaultOptions(sys_M_RF4);
        sys_M_RF16 = SystemSettings(M, N_RF_16, SNR_default);
        options_M_RF16 = MyDefaultOptions(sys_M_RF16);
        %% display the simulation progress
        disp(['M=', num2str(M), '  trial=', num2str(trial_idx)]);

        %% system and algorithm parameters

        y_M_RF4 = y_RF4(1:M);
        y_nf_RF4 = Y_nf_RF4(1:M);
        y_M_RF16 = y_RF16(1:M);
        y_nf_RF16 = Y_nf_RF16(1:M);
        W_M = W(:,1:M);
        n_M = n(1:M);

        %% conventional CRAF with no joint phase
        [h_CRAF] = CRAF_CE(y_noncoherent(1:M), W_Gaussian(:,1:M), sys_M_RF4, options_M_RF4, h);
        SE = SE_rotate(h_CRAF, h);
        SE_CRAF(M_idx, trial_idx) = SE;
        SNRLoss_CRAF(M_idx, trial_idx) = abs((sign(h_CRAF)'*h)/(sign(h)'*h));

        %% partially coherent compressive channel estimation
        %% N_RF = 4
        h0 = InitGuess_UPA(y_M_RF4, W_M, sys_M_RF4, options_M_RF4);
        [h_PCCPR_list] = PCCPR_UPA(y_M_RF4, W_M, h0, sys_M_RF4, options_M_RF4);
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_RF4(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_RF4(M_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_M_RF4.T))'*h)/(sign(h)'*h));

        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y_M_RF4, W_M, h0, sys_M_RF4, options_M_RF4);
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_RF4(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_RF4(M_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_M_RF4.T2))'*h)/(sign(h)'*h));

        %% N_RF = 16
        h0 = InitGuess_UPA(y_M_RF16, W_M, sys_M_RF16, options_M_RF16, sys_M_RF16.N*noise_var);
        [h_PCCPR_list] = PCCPR_UPA(y_M_RF16, W_M, h0, sys_M_RF16, options_M_RF16);
        SE_list = SE_rotate(h_PCCPR_list, h);
        SE_PCCPR_RF16(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_RF16(M_idx, trial_idx) = abs((sign(h_PCCPR_list(:,sys_M_RF16.T))'*h)/(sign(h)'*h));

        [h_PCCPR_SR_list] = PCCPR_SR_UPA(y_M_RF16, W_M, h0, sys_M_RF16, options_M_RF16);
        SE_list = SE_rotate(h_PCCPR_SR_list, h);
        SE_PCCPR_SR_RF16(M_idx, trial_idx, :) = SE_list;
        SNRLoss_PCCPR_SR_RF16(M_idx, trial_idx) = abs((sign(h_PCCPR_SR_list(:,sys_M_RF16.T2))'*h)/(sign(h)'*h));

        %% NOMP with perfect phase information
        y_M = W_M'*h + n_M;
        h_est_NOMP = NOMP_2D('2nd', sys_M_RF4.N1, sys_M_RF4.N2, W_M, y_M, options_M_RF4.Np, 1e-5, 5);
        SE_NOMP(M_idx, trial_idx) = norm(h_est_NOMP - h, 2)^2;
        SNRLoss_NOMP(M_idx, trial_idx) = abs((sign(h_est_NOMP)'*h)/(sign(h)'*h));

    end
end

for M_idx = 1:M_list_size
    %% calculate the success rates
    success_rate_CRAF(M_idx) = mean(SNRLoss_CRAF(M_idx,:)>0.5);
    success_rate_PCCPR_RF4(M_idx) = mean(SNRLoss_PCCPR_RF4(M_idx,:)>0.5);
    success_rate_PCCPR_SR_RF4(M_idx) = mean(SNRLoss_PCCPR_SR_RF4(M_idx,:)>0.5);
    success_rate_PCCPR_RF16(M_idx) = mean(SNRLoss_PCCPR_RF16(M_idx,:)>0.5);
    success_rate_PCCPR_SR_RF16(M_idx) = mean(SNRLoss_PCCPR_SR_RF16(M_idx,:)>0.5);
    success_rate_NOMP(M_idx) = mean(SNRLoss_NOMP(M_idx,:)>0.5);
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
plot(M_list,success_rate_PCCPR_RF4,'b-o','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_RF4,'c-*','LineWidth',1.5);
hold on;
plot(M_list,success_rate_NOMP(:),'g-v','LineWidth',1.5);
hold on;
plot(M_list,success_rate_CRAF(:),'r-^','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_RF16,'b-o','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_RF16,'c-*','LineWidth',1.5);
legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF');
grid on;
xlabel('M');
ylabel('Channel estimation success rate');

figure;
plot(M_list,mean(SNRLoss_PCCPR_RF4, 2),'b-o','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_SR_RF4, 2),'c-*','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_NOMP, 2),'g-v','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_CRAF, 2),'r-^','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_RF16, 2),'b-o','LineWidth',1.5);
hold on;
plot(M_list,mean(SNRLoss_PCCPR_SR_RF16, 2),'c-X','LineWidth',1.5);
legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF');
grid on;
xlabel('M');
ylabel('mean SNR loss');

figure;
semilogy(M_list,NMSE_PCCPR_1,'b-o','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_1,'c-*','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_NOMP,'g-v','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_CRAF,'r-^','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_2,'b-o','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_2,'c-*','LineWidth',1.5);
legend('PCCPR','PCCPR-SR','NOMP with perfectly known phase','CRAF');
grid on;
xlabel('M');
ylabel('NMSE performance');