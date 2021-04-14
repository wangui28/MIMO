clc;clear all;close all;

sys=SystemSet();
Nt = sys.Nt; %发射天线数
Nr = sys.Nr; %接收天线数
Nx = sys.Nx; %发射导频数
Ny = sys.Ny; %接收导频数
L = sys.L; %路径数

N_trial = 2;

SNR_list = -30:5:0; %信噪比取值
SNR_list_size = numel(SNR_list);

h_norm = zeros(N_trial, 1);
NMSE_PCCPR_SR = zeros(SNR_list_size, N_trial, sys.T2+1, 1);

%% try Monte-Carlo simulations for N_trial instances
for trial_idx = 1:N_trial
    
	%% 导频编码矩阵
	X = exp(1i*2*pi*rand(Nt, Nx)); 
        
	%% 组合解码矩阵
	W = exp(1i*2*pi*rand(Nr, Ny));
        
	C = kron(X.',W'); %NxNy * NtNr
    
%     % generate random structured phase offsets
%     phase_1 = 2*pi*rand(1);
%     phase_11 = zeros(sys.M,1);
%     for i=1:sys.M
%         phase_11(i)=phase_1;
%     end
%     phase = diag(phase_11);
    phase_l = 2*pi*rand(sys.block, 1);
    phase = zeros(sys.M, 1);
    for block_idx = 1:sys.block
        phase((block_idx-1)*sys.size+1:block_idx*sys.size,1) = phase_l(block_idx, 1)*ones(sys.size, 1);
    end
    
    [h, alpha, theta_t, theta_r] = Channel(Nr, Nt, L); % 信道
    h_norm(trial_idx, :) = norm(h);
    
    %% 对每一个信噪比值循环
    for SNR_idx = 1:SNR_list_size
        
        SNR = SNR_list(SNR_idx);
    
        disp(['trial=', num2str(trial_idx), '    SNR=', num2str(SNR)]);  %窗口显示仿真进度
    
        %% 初始化噪声
        noise = 10^(-SNR/10)/Nx/Ny;
        n = sqrt(noise/2)*(normrnd(0, 1, Ny*Nx, Ny*Nx) + 1i*normrnd(0, 1, Ny*Nx, Ny*Nx)); %复噪声矩阵
        
        y = (C*h).*exp(1i*phase(:,1)) + n(:,1);
        
        h0 = Initialization(y, C, sys);
        [h_PCCPR_SR_list] = PCCPR_SR_ULA(y, C, h0, sys);
        
        NMSE_list = NMSE_rotate(h_PCCPR_SR_list, h);
        NMSE_PCCPR_SR(SNR_idx, trial_idx, :, :) = NMSE_list;
    end     
end

NMSE_PCCPR_SR = sum(sum(squeeze(NMSE_PCCPR_SR(:,:,end,:)),2),3)/sum(sum(h_norm.^2));

figure;
semilogy(SNR_list,NMSE_PCCPR_SR,'b-*','LineWidth',1.5);
grid on;
xlabel('SNR(dB)'); ylabel('NMSE');