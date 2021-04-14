close all;clear all;clc;

Nt0 = [8;8];
Nt = Nt0(1)*Nt0(2);  
Nr0 = [8;8];
Nr = Nt0(1)*Nt0(2);
Ny = 32;
Nx = 32;
N_RF = 4;
L = 3;
d = 0.5;
H = zeros(Nr, Nt);

X = 1/sqrt(Nt)*exp(-1i*2*pi*rand(Nt, Nx));
mu_X = max(max(abs(X*X'-Nx/Nt*eye(Nt))));
W = 1/sqrt(Nr)*exp(-1i*2*pi*rand(Ny, Nr));
mu_W = max(max(abs(W'*W-Ny/Nr*eye(Nr))));
for i = 1:1000000
    X0 = 1/sqrt(Nt)*exp(-1i*2*pi*rand(Nt, Nx));
    mu_X0 = max(max(abs(X0*X0'-Nx/Nt*eye(Nt))));
    if (mu_X0 < mu_X)
        X = X0;
        mu_X = mu_X0;
    end
    W0 = 1/sqrt(Nr)*exp(-1i*2*pi*rand(Ny, Nr));
    mu_W0 = max(max(abs(W0'*W0-Ny/Nr*eye(Nr))));
    if (mu_W0 < mu_W)
        W = W0;
        mu_W = mu_W0;
    end
end

SNR_list = -5:2.5:20;
sample_num = 500;
nmse_result = zeros(numel(SNR_list),sample_num);
alpha_true = zeros(numel(SNR_list),sample_num,L);
phi_t_true = zeros(numel(SNR_list),sample_num,L,2);
phi_r_true = zeros(numel(SNR_list),sample_num,L,2);
L_result = zeros(numel(SNR_list),sample_num);
theta_result = zeros(numel(SNR_list),sample_num,4,10);
z_result = zeros(numel(SNR_list),sample_num,10);

for sample_ii = 1:sample_num
    H = zeros(Nr, Nt);
    alpha = zeros(L,1);
    alpha(1:L) = (normrnd(0, 1, L, 1) + 1i*normrnd(0, 1, L, 1)) / sqrt(2);
    while (find(abs(alpha)<0.01))
        alpha(1:L) = (normrnd(0, 1, L, 1) + 1i*normrnd(0, 1, L, 1)) / sqrt(2);
    end
    alpha = sort(alpha, 'descend');
    phi_t = 2*rand(L,2)-1;%virtual AoD
    phi_r = 2*rand(L,2)-1;%virtual AoA
        
    for l = 1:L
        at = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*d*phi_t(l,1)), exp(-1i*2*pi*[0:Nt0(2)-1]'*d*phi_t(l,2)));
        ar = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*d*phi_r(l,1)), exp(-1i*2*pi*[0:Nr0(2)-1]'*d*phi_r(l,2)));
        H = H + alpha(l)*ar*at';
    end
    
    for snr_ii = 1:numel(SNR_list)
        snr = SNR_list(snr_ii);
        noise = sqrt(10^(-snr/10)/2);
        
        H = zeros(Nr, Nt);
        alpha = zeros(L,1);
        alpha(1:L) = (normrnd(0, 1, L, 1) + 1i*normrnd(0, 1, L, 1)) / sqrt(2);
        while (find(abs(alpha)<0.01))
            alpha(1:L) = (normrnd(0, 1, L, 1) + 1i*normrnd(0, 1, L, 1)) / sqrt(2);
        end
        alpha = sort(alpha, 'descend');
        phi_t = 2*rand(L,2)-1;%virtual AoD
        phi_r = 2*rand(L,2)-1;%virtual AoA
        
        for l = 1:L
            at = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*d*phi_t(l,1)), exp(-1i*2*pi*[0:Nt0(2)-1]'*d*phi_t(l,2)));
            ar = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*d*phi_r(l,1)), exp(-1i*2*pi*[0:Nr0(2)-1]'*d*phi_r(l,2)));
            H = H + alpha(l)*ar*at';
        end

        Y = W*(H*X + noise*(normrnd(0, 1, Nr, Nx) + 1i*normrnd(0, 1, Nr, Nx)));
        Rth = noise*sqrt(Ny*Nx);
        [theta_es,z_es,err]=IR_SURE_CE_UPA(Y,X,W,Nx,Nt0,Nr0,Ny,Rth);
        H_es = zeros(Nr, Nt);
        at = zeros(Nt,1);
        ar = zeros(Nr,1);
        for l = 1:numel(z_es)
            at = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta_es(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta_es(2,l)));
            ar = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta_es(3,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta_es(4,l)));
            H_es = H_es + z_es(l)*ar*at';
        end
        
        nmse_sample = sum(sum(abs(H-H_es).^2))/sum(sum(abs(H).^2));
        disp(['sample_ii=' num2str(sample_ii) ' snr=' num2str(SNR_list(snr_ii)) ' nmse=' num2str(nmse_sample) ' err=' num2str(err)]);
        
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
        phi_t_true(snr_ii, sample_ii,:,:) = phi_t;
        phi_r_true(snr_ii, sample_ii,:,:) = phi_r;
        nmse_result(snr_ii, sample_ii) = nmse_sample;
        L_sample = numel(z_es);
        L_result(snr_ii, sample_ii) = L_sample;
        theta_result(snr_ii, sample_ii, :, 1:L_sample) = theta_es(:,:);
        z_result(snr_ii, sample_ii, 1:L_sample) = z_es(:);
    end
    save result.mat
end