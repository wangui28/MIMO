close all;clear all;clc;
timer = 0;
Nt = 64;
Nr = 64;
Ny = 36;
Nx = 36;
N_RF = 4;
L = 3;
d = 0.5;
H = zeros(Nr, Nt);

SNR_list = -5:2.5:0;
sample_num = 100;
nmse_result = zeros(numel(SNR_list),sample_num);
alpha_true = zeros(numel(SNR_list),sample_num,L);
phi_t_true = zeros(numel(SNR_list),sample_num,L);
phi_r_true = zeros(numel(SNR_list),sample_num,L);
L_result = zeros(numel(SNR_list),sample_num);
theta_result = zeros(numel(SNR_list),sample_num,2,10);
z_result = zeros(numel(SNR_list),sample_num,10);
for snr_ii = 1:numel(SNR_list)
    snr = SNR_list(snr_ii);
    noise = sqrt(10^(-snr/10)/2);
    for sample_ii = 1:sample_num
        H = zeros(Nr, Nt);
        alpha = zeros(L,1);
        alpha(1) = exp(1i*2*pi*rand(1));
        alpha(2:L) = (normrnd(0, 0.1, L-1, 1) + 1i*normrnd(0, 0.1, L-1, 1)) / sqrt(2);
        while (find(abs(alpha)<0.01))
            alpha(2:L) = (normrnd(0, 0.1, L-1, 1) + 1i*normrnd(0, 0.1, L-1, 1)) / sqrt(2);
        end
        alpha = sort(alpha, 'descend');
        phi_t = 2*rand(L,1)-1;%virtual AoD
        phi_r = 2*rand(L,1)-1;%virtual AoA
        
        for l = 1:L
            at = exp(-1i*2*pi*[0:Nt-1]'*d*phi_t(l));
            ar = exp(-1i*2*pi*[0:Nr-1]'*d*phi_r(l));
            H = H + alpha(l)*ar*at';
        end

        X = 1/sqrt(Nt)*exp(-1i*2*pi*rand(Nt, Nx));
        W = 1/sqrt(Nr)*exp(-1i*2*pi*rand(Ny, Nr));
        Y = W*(H*X + noise*(normrnd(0, 1, Nr, Nx) + 1i*normrnd(0, 1, Nr, Nx)));
        Rth = noise*sqrt(Ny*Nx);
        tic;
        [theta_es,z_es,err]=IR_SURE_CE(Y,X,W,Nx,Nt,Nr,Ny,Rth);
        timer = timer + toc;
        H_es = zeros(Nr, Nt);
        at = zeros(Nt,1);
        ar = zeros(Nr,1);
        for l = 1:numel(z_es)
            at = exp(-1i*2*pi*[0:Nt-1]'*theta_es(1,l));
            ar = exp(-1i*2*pi*[0:Nr-1]'*theta_es(2,l));
            H_es = H_es + z_es(l)*ar*at';
        end
        
        nmse_sample = sum(sum(abs(H-H_es).^2))/sum(sum(abs(H).^2));
        disp(['snr=' num2str(SNR_list(snr_ii)) ' sample_ii=' num2str(sample_ii) ' nmse=' num2str(nmse_sample) ' err=' num2str(err)]);
        
       %% spectral efficiency
        [U_perfectCSI,S,V_perfectCSI] = svd(H);
        P_perfectCSI = V_perfectCSI(:,1:N_RF);
        Q_perfectCSI = U_perfectCSI(:,1:N_RF);
        R_perfectCSI = 0.5*noise*noise*(Q_perfectCSI'*Q_perfectCSI);
        SE_perfectCSI = log2(det(eye(N_RF)+(1/N_RF)*inv(R_perfectCSI)*((Q_perfectCSI'*H*P_perfectCSI)*(Q_perfectCSI'*H*P_perfectCSI)')));
        
        [U_estCSI_SURE,S,V_estCSI_SURE] = svd(H_es);
        P_estCSI_SURE = V_estCSI_SURE(:,1:N_RF);
        Q_estCSI_SURE = U_estCSI_SURE(:,1:N_RF);
        R_estCSI_SURE = 0.5*noise*noise*(Q_estCSI_SURE'*Q_estCSI_SURE);
        SE_estCSI_SURE = log2(det(eye(N_RF)+(1/N_RF)*inv(R_estCSI_SURE)*((Q_perfectCSI'*H*P_estCSI_SURE)*(Q_perfectCSI'*H*P_estCSI_SURE)')));
        
        alpha_true(snr_ii, sample_ii,:) = alpha;
        phi_t_true(snr_ii, sample_ii,:) = phi_t;
        phi_r_true(snr_ii, sample_ii,:) = phi_r;
        nmse_result(snr_ii, sample_ii) = nmse_sample;
        L_sample = numel(z_es);
        L_result(snr_ii, sample_ii) = L_sample;
        theta_result(snr_ii, sample_ii, :, 1:L_sample) = theta_es(:,:);
        z_result(snr_ii, sample_ii, 1:L_sample) = z_es(:);
    end
    save result.mat
end