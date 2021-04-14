clc; clear all; close all

N1=16; N2=16; N=N1*N2; % the number of antennas
M=128; % the number of measurements
L=3; % the number of paths
d=0.5; % the space of antennas
type=1; % ULA or UPA
K=2000;  % the number of samples

% DFT matrix
if type==1
    UN=(1/sqrt(N))*exp(1i*2*pi*[-(N-1):2:N-1]'/2*d*[-(N-1):2:N-1]*(1/N)); % the DFT matrix
end
if type==2
    UN1=(1/sqrt(N1))*exp(1i*2*pi*[-(N1-1):2:N1-1]'/2*d*[-(N1-1):2:N1-1]*(1/N1));
    UN2=(1/sqrt(N2))*exp(1i*2*pi*[-(N2-1):2:N2-1]'/2*d*[-(N2-1):2:N2-1]*(1/N2));
    UN=kron(UN1,UN2); % the DFT matrix
end

x=zeros(N,K); % the beamspace channel 
for k=1:K
    x(:,k)=UN.'*generate_channel(N1,N2,L,type);
end

load(['CSmatrix256',num2str(M),'.mat'])

SNR_dB=[0:5:20];
SNR_linear=10.^(SNR_dB/10.);
OMP_nmse=zeros(1,length(SNR_dB));
AMP_nmse=zeros(1,length(SNR_dB));
LAMP_nmse=zeros(1,length(SNR_dB));
GM_LAMP_nmse=zeros(1,length(SNR_dB));

for iS=1:length(SNR_dB)
    sigma2=1/SNR_linear(iS);   
    noise = sqrt(sigma2)*(randn(N,K)+1i*randn(N,K))/sqrt(2);
    y=A*(x+noise);
    % OMP algorithm 
    OMP_xhat=[];
    for i=1:K
        yi=y(:,i);
        OMP_xihat=OMP(yi,A,24);
        OMP_xhat=[OMP_xhat,OMP_xihat];
    end
   OMP_nmse(iS)=10*log10(mean((sum(abs(OMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
   % AMP algorithm
   [AMP_xhat]=AMP(y,A);
   AMP_nmse(iS)=10*log10(mean((sum(abs(AMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
   
   if SNR_dB(iS)<10
       snr_range=1;
   end
   if SNR_dB(iS)>=10
       snr_range=2;
   end
   % LAMP algorithm
   [LAMP_xhat]=LAMP(y,A,type,snr_range);
   LAMP_nmse(iS)=10*log10(mean((sum(abs(LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
   % GM-LAMP algorithm
   [GM_LAMP_xhat]=GM_LAMP(y,A,type,snr_range);
   GM_LAMP_nmse(iS)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))

end



figure('color',[1,1,1])

plot(SNR_dB,OMP_nmse,'ms-',...
    SNR_dB,AMP_nmse,'g*-',...
    SNR_dB,LAMP_nmse,'bd-',...  
    SNR_dB,GM_LAMP_nmse,'r^-','LineWidth',1.5);

L1=legend('OMP [11]','AMP [15]','LAMP [29]','Proposed GM-LAMP','Location','NorthEast');set(gca,'FontSize',11, 'FontName','Arial')
set(L1,'FontSize',11, 'FontName','Arial')
xlabel('SNR (dB)','FontSize',11,'FontName','Arial')
ylabel('NMSE (dB)','FontSize',11,'FontName','Arial')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
% set(gca,'XTick',[0:5:20])
% set(gca,'YTick',[-21:4:1])
% axis([0,20,-21,1])
