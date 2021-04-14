clear all;close all;clc;
load NRF_result.mat


figure;
plot(N_RF_list,mean(-10*log10(SNRLoss_PCCPR_1), 2),'b-^','LineWidth',1.5);
hold on;
plot(N_RF_list,mean(-10*log10(SNRLoss_PCCPR_SR_1), 2),'c-v','LineWidth',1.5);
hold on;
plot(N_RF_list,mean(-10*log10(SNRLoss_CRAF), 2),'r-x','LineWidth',1.5,'Markersize',10);
hold on;
plot(N_RF_list,mean(-10*log10(SNRLoss_OMP), 2),'m-*','LineWidth',1.5,'Markersize',10);
hold on;
plot(N_RF_list,mean(-10*log10(SNRLoss_NOMP), 2),'g-d','LineWidth',1.5);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','CRAF [20]','OMP assuming known phase [8]','NOMP assuming known phase [34]')
grid on;
xlim([0 18]);
xlabel('{\it N}_{\rm RF}');
ylabel('Mean SNR loss (dB)');
set(gcf,'position',[200,200,520,390]);

figure;
semilogy(N_RF_list,NMSE_PCCPR_1,'b-^','LineWidth',1.5);
hold on;
semilogy(N_RF_list,NMSE_PCCPR_SR_1,'c-v','LineWidth',1.5);
hold on;
semilogy(N_RF_list,NMSE_CRAF,'r-x','LineWidth',1.5,'Markersize',10);
hold on;
semilogy(N_RF_list,NMSE_OMP,'m-*','LineWidth',1.5,'Markersize',10);
hold on;
semilogy(N_RF_list,NMSE_NOMP,'g-d','LineWidth',1.5);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','CRAF [20]','OMP assuming known phase [8]','NOMP assuming known phase [34]')
grid on;
xlim([0 18]);
ylim([0.001 1]);
xlabel('{\it N}_{\rm RF}');
ylabel('NMSE performance');
set(gcf,'position',[200,200,520,390]);



clear all;
load init_result.mat
figure;
plot(M_list,success_rate_PCCPR_opt,'b-^','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_opt,'c-v','LineWidth',1.5);
hold on;
p=plot(M_list,success_rate_PCCPR_bm,'g-x','LineWidth',1.5,'markersize',10);
set(p,'Color',[0 0.6 0]);
hold on;
plot(M_list,success_rate_PCCPR_SR_bm,'g-*','LineWidth',1.5,'markersize',10);
hold on;
p=plot(M_list,success_rate_PCCPR_rnd,'r-s','LineWidth',1.5);
set(p,'Color',[0.6 0 1]);
hold on;
plot(M_list,success_rate_PCCPR_SR_rnd,'r-d','LineWidth',1.5);
legend('On-grid PC-CPR, proposed initialization','Off-grid PC-CPR, proposed initialization','On-grid PC-CPR, the initialization in [20]','Off-grid PC-CPR, the initialization in [20]','On-grid PC-CPR, random initialization','On-grid PC-CPR, random initialization');
grid on;
xlabel('{\it M}');
ylabel('Channel estimation success rate');
set(gcf,'position',[200,200,520,390]);






clear all;
load M_result.mat

figure;
plot(M_list,success_rate_PCCPR_RF4,'b-^','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_RF4,'c-v','LineWidth',1.5);
hold on;
plot(M_list,success_rate_CRAF(:),'r-x','LineWidth',1.5,'markersize',10);
hold on;
plot(M_list,success_rate_NOMP(:),'g-d','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_RF16,'b-^','LineWidth',1.5);
hold on;
plot(M_list,success_rate_PCCPR_SR_RF16,'c-v','LineWidth',1.5);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','CRAF [20]','NOMP assuming known phase [34]');
grid on;
xlabel('{\it M}');
ylabel('Channel estimation success rate');
set(gcf,'position',[200,200,520,390]);


figure;
semilogy(M_list,NMSE_PCCPR_1,'b-^','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_1,'c-v','LineWidth',1.5);
hold on;
semilogy(M_list,min(NMSE_CRAF,1),'r-x','LineWidth',1.5,'markersize',10);
hold on;
semilogy(M_list,NMSE_NOMP,'g-d','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_2,'b-^','LineWidth',1.5);
hold on;
semilogy(M_list,NMSE_PCCPR_SR_2,'c-v','LineWidth',1.5);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','CRAF [20]','NOMP assuming known phase [34]');
grid on;
ylim([0 10]);
xlabel('{\it M}');
ylabel('NMSE performance');
set(gcf,'position',[200,200,520,390]);



clear all;
load SNR_result.mat

figure;
semilogy(SNR_list,NMSE_PCCPR_1,'b-^','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_PCCPR_SR_1,'c-v','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_CRAF,'r-x','LineWidth',1.5,'Markersize',10);
hold on;
semilogy(SNR_list,NMSE_NOMP,'g-d','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_PCCPR_2,'b-^','LineWidth',1.5);
hold on;
semilogy(SNR_list,NMSE_PCCPR_SR_2,'c-v','LineWidth',1.5);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','CRAF [20]','NOMP assuming known phase [34]')
grid on;
xlabel('SNR (dB)');
ylabel('NMSE performance');
set(gcf,'position',[200,200,520,390]);






clear all;
load MU_result.mat
figure;
plot(SNR_list,mean(SumRate_PCCPR,2),'b-^','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SumRate_PCCPR_SR,2),'c-v','LineWidth',1.5);
hold on;
plot(SNR_list,mean(SumRate_CRAF,2),'r-x','LineWidth',1.5,'markersize',10);
hold on;
plot(SNR_list,mean(SumRate_NOMP,2),'g-d','LineWidth',1.5);
set(gcf,'position',[200,200,520,390]);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','CRAF [20]','NOMP assuming known phase [34]');
grid on;
xlabel('P_T/\sigma_n^2 (dB)');
ylabel('Multiuser sum-rate (bit/s/Hz)');



clear all;
load noncoherentCE_SNR0PATH3Rician.mat
load SNR0PATH3Rician.mat

r = 3;
success_th = 0.5;

figure;
plot(Mlist,mean(SNRLoss_PCCPR>success_th,2),'b-^','LineWidth',1.5);
hold on;
plot(Mlist,mean(SNRLoss_PCCPR_SR>success_th,2),'c-v','LineWidth',1.5);
hold on;
p=plot(Mlist,mean(SNRLoss_FBA>success_th,2),'c-o','LineWidth',1.5);
set(p,'Color',[0.6 0 1]);
hold on;
p=plot(Mlist,mean(SNRLoss(r:r:r*Mlist_size,:)>0.35,2),'g-s','LineWidth',1.5);
set(p,'Color',[0 0.6 0]);
hold on;
plot(Mlist,mean(SNRLoss_NOMP>success_th,2),'g-d','LineWidth',1.5);
legend('Algorithm 1 (on-grid PC-CPR)','Algorithm 2 (off-grid PC-CPR)','FBA [14]','Noncoherent CE [17]','NOMP assuming known phase [34]')
grid on;
xlabel('{\it M}');
ylabel('Channel estimation success rate');
set(gcf,'position',[200,200,520,390]);