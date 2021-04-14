clc; clear all; close all

N=256; % the number of antennas
M=128; % the number of measurements
L=3; % the number of paths 
d=0.5; % the space of antennas
K=2000;  % the number of samples

UN=(1/sqrt(N))*exp(1i*2*pi*[-(N-1):2:N-1]'/2*d*[-(N-1):2:N-1]*(1/N)); % the DFT matrix

x=zeros(N,K); % the beamspace channel 
for k=1:K
    x(:,k)=UN.'*generate_channel(16,16,L,1);
end

load(['CSmatrix256',num2str(M),'.mat'])

SNR_dB=[0:2:20];
SNR_linear=10.^(SNR_dB/10.);
nmse_5dB=zeros(1,length(SNR_dB));
nmse_10dB=zeros(1,length(SNR_dB));
nmse_15dB=zeros(1,length(SNR_dB));
nmse_0to10dB=zeros(1,length(SNR_dB));
nmse_10to20dB=zeros(1,length(SNR_dB));
nmse_0to20dB=zeros(1,length(SNR_dB));


for iS=1:length(SNR_dB)
    sigma2=1/SNR_linear(iS);   
    noise = sqrt(sigma2)*(randn(N,K)+1i*randn(N,K))/sqrt(2);
    y=A*(x+noise);

   % for the trained GM-LAMP network at 5 dB
   trainedfile_name='Trained_SV_ULA_gm_2561285dB.mat';
   [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name);
   nmse_5dB(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
      
   % for the trained GM-LAMP network at 10 dB
   trainedfile_name='Trained_SV_ULA_gm_25612810dB.mat';
   [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name);
   nmse_10dB(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
      
   % for the trained GM-LAMP network at 15 dB
   trainedfile_name='Trained_SV_ULA_gm_25612815dB.mat';
   [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name);
   nmse_15dB(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
   
   % for the trained GM-LAMP network from 0 dB to 10 dB
   trainedfile_name='Trained_SV_ULA_gm_2561280to10dB.mat';
   [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name);
   nmse_0to10dB(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
      
   % for the trained GM-LAMP network from 0 dB to 20 dB
   trainedfile_name='Trained_SV_ULA_gm_2561280to20dB.mat';
   [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name);
   nmse_0to20dB(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
      
   % for the trained GM-LAMP network form 10 dB to 20 dB
   trainedfile_name='Trained_SV_ULA_gm_25612810to20dB.mat';
   [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name);
   nmse_10to20dB(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))

end


figure('color',[1,1,1])

plot(SNR_dB,nmse_5dB,'bs--',...  
    SNR_dB,nmse_10dB,'b^--',...
    SNR_dB,nmse_15dB,'bo--',...
    SNR_dB,nmse_0to10dB,'rs-',...
    SNR_dB,nmse_0to20dB,'r^-',...
    SNR_dB,nmse_10to20dB,'ro-','LineWidth',1.5);

L1=legend('SSNR:  5 dB','SSNR: 10 dB','SSNR: 15 dB','MSNR:  0-10 dB','MSNR:  0-20 dB','MSNR: 10-20 dB','Location','NorthEast');set(gca,'FontSize',11, 'FontName','Arial')
set(L1,'FontSize',11, 'FontName','Arial')
xlabel('SNR (dB)','FontSize',11,'FontName','Arial')
ylabel('NMSE (dB)','FontSize',11,'FontName','Arial')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'XTick',[0:2:20])
set(gca,'YTick',[-23:5:5])
axis([0,20,-23,5])