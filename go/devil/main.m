clc;close all;clear all;

Nt = 32; % 发射天线数
Nr = 32; % 接收天线数
Nx = 16; % 发射导频数
Ny = 16; % 接收导频数
L = 5; % 路径数
sys = SystemSet(Nt, Nr, Nx, Ny, L);
alg = AlgorithmSet(sys);

SNR_list = -30:5:0; % SNR
SNR_list_size = numel(SNR_list);

N_trial = 2;
K = sys.Nx;

%% the containers to put the results
h_norm = zeros(N_trial, K);

NMSE_PCCPR_SR = zeros(SNR_list_size, N_trial, alg.T_off+1, K);

%% Monte-Carlo Simulations for N_trial instances
for trial_idx = 1:N_trial
         
    C = exp(1i*2*pi*rand(sys.N, sys.M));
    
    % 相位偏移
    phase_l = 2*pi*rand(sys.Nx, K);
    phase = zeros(sys.M, K);
    for k_idx = 1:K
        for block_idx = 1:sys.Nx
            phase((block_idx-1)*sys.Ny+1:block_idx*sys.Ny,k_idx) = phase_l(block_idx, k_idx)*ones(sys.Ny, 1);
        end
    end
%     phase_noncoherent = 2*pi*rand(sys.M, K);
    
    % 信道建模
    H = zeros(sys.N, K);
    for k_idx = 1:K
        [h, alpha, theta_t, theta_r] = Channel(sys.Nr, sys.Nt, sys.L);
        H(:,k_idx) = h;
        h_norm(trial_idx, k_idx) = norm(h);
    end
    
    %% go through all SNR's
    for SNR_idx = 1:SNR_list_size
        
        SNR = SNR_list(SNR_idx);
    
        disp(['trial=', num2str(trial_idx), '    SNR=', num2str(SNR)]);  %窗口显示仿真进度
    
        % 噪声
        noise = 10^(-SNR/10)/sys.Nx;
        n = sqrt(noise/2)*(normrnd(0, 1, sys.M, K)+1i*normrnd(0, 1, sys.M, K)); %复噪声矩阵
        
        for k_idx = 1:K
            h = H(:,k_idx);
            y = (C'*h).*exp(1i*phase(:,k_idx)) + n(:,k_idx);
        
            h0 = Initialization(y, C, sys, alg);
            [h_PCCPR_SR_list] = PCCPR_OffGrid(y, C, h0, sys, alg);
        
            NMSE_list = NMSE_rotate(h_PCCPR_SR_list, h);
            NMSE_PCCPR_SR(SNR_idx, trial_idx, :, k_idx) = NMSE_list;
        end
    end     
end

%% calculate the NMSE
NMSE_PCCPR_SR = sum(sum(squeeze(NMSE_PCCPR_SR(:,:,end,:)),2),3)/sum(sum(h_norm.^2));

figure;
semilogy(SNR_list,NMSE_PCCPR_SR,'b-*','LineWidth',1.5);
xlabel('SNR(dB)'); ylabel('NMSE');grid on;