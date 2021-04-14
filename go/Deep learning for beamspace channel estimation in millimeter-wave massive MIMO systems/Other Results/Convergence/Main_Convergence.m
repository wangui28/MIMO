clc; clear all; close all

N=256; % the number of antennas
M=128; % the number of measurements
L=3; % the number of paths 
d=0.5; % the space of antennas
K=2000;  % the number of samples
T_all=1:12; % the number of layers

UN=(1/sqrt(N))*exp(1i*2*pi*[-(N-1):2:N-1]'/2*d*[-(N-1):2:N-1]*(1/N)); % the DFT matrix

x=zeros(N,K); % the beamspace channel 
for k=1:K
    x(:,k)=UN.'*generate_channel(16,16,L,1);
end

load(['CSmatrix256',num2str(M),'.mat'])


T_5dB=zeros(1,length(T_all));
T_15dB=zeros(1,length(T_all));

for iT=1:length(T_all)
    SNR_linear=10.^(5/10.);
    sigma2=1/SNR_linear;   
    noise = sqrt(sigma2)*(randn(N,K)+1i*randn(N,K))/sqrt(2);
    y=A*(x+noise);
    T=T_all(iT);
    trainedfile_name=['Trained_SV_ULA_gm_256128_for_',num2str(5),'dB_T=',num2str(T),'.mat'];
    [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name,T);
    T_5dB(iT)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
    
    SNR_linear=10.^(15/10.);
    sigma2=1/SNR_linear;   
    noise = sqrt(sigma2)*(randn(N,K)+1i*randn(N,K))/sqrt(2);
    y=A*(x+noise);
    T=T_all(iT);
    trainedfile_name=['Trained_SV_ULA_gm_256128_for_',num2str(15),'dB_T=',num2str(T),'.mat'];
    [GM_LAMP_xhat]=GM_LAMP(y,A,trainedfile_name,T);
    T_15dB(iT)=10*log10(mean((sum(abs(GM_LAMP_xhat-x).^2,1)))./mean(sum(abs(x).^2,1)))
end


figure('color',[1,1,1]) 
plot(T_all,T_5dB,'rd-',...
    T_all,T_15dB,'r^-','LineWidth',1.5);
L1=legend(' 5 dB','15 dB','Location','NorthEast');

set(gca,'FontSize',11, 'FontName','Arial')
set(L1,'FontSize',11, 'FontName','Arial')
xlabel('Layer \it T','FontSize',11,'FontName','Arial')
ylabel('NMSE (dB)','FontSize',11,'FontName','Arial')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'XTick',[1:1:12])
set(gca,'YTick',[-23:3:-8])
axis([1,12,-23,-8])