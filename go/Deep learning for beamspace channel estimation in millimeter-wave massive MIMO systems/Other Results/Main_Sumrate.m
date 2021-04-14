clc; clear all; close all

NMSE=[0:-5:-30];  
NMSE_linear=10.^(NMSE/10.);
N_iter=2000; 
N = 256; % number of beams (transmit antennas)
K = 16; % number of users
L = 3; % number of paths per user
lamada = 1; % wavelength
d = lamada/2; % the space of antennas
down_SNR=10.^(10/10.);  % 10dB


UN=(1/sqrt(N))*exp(1i*2*pi*[-(N-1):2:N-1]'/2*d*[-(N-1):2:N-1]*(1/N));

for i_e=1:length(NMSE) 
    i_e
    sigma2=NMSE_linear(i_e);
    temp = 0; temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0; temp5 = 0; temp6 = 0;
    error1 = 0; error2 = 0; error3 = 0; error4 = 0; 
    for iter = 1:N_iter
        H=zeros(N,K);
        for k=1:K
            H(:,k)= generate_channel(16,16,L,1); % generate the signalspace channel
        end
        H_beam = UN.'*H; % beamspace channel
        
        %%%%%%% perfect CSI, IA-BS
        [Hr,Hr_e] = IA_BS(H_beam,H_beam); 
        F = Hr_e*inv(Hr_e'*Hr_e);
        beta = sqrt(K/trace(F*F'));
        H_eq=Hr'*F;
        for k=1:K
            sum_inf=sum(abs(H_eq(:,k)).^2)-abs(H_eq(k,k))^2;
            temp=temp+log2(1+abs(H_eq(k,k))^2/(sum_inf+K/(down_SNR*beta^2)));
        end
       
        %%%%%% nmse:-10dB CSI, IA-BS
        noise1= sqrt(sigma2)*(randn(N,K)+1i*randn(N,K))/sqrt(2);
        H_beam1=H_beam+noise1;
        [Hr1,Hr1_e] = IA_BS(H_beam1,H_beam); 
        F1 = Hr1_e*inv(Hr1_e'*Hr1_e);
        beta1 = sqrt(K/trace(F1*F1'));
        H_eq1=Hr1'*F1;
        for k=1:K
            sum_inf=sum(abs(H_eq1(:,k)).^2)-abs(H_eq1(k,k))^2;
            temp1=temp1+log2(1+abs(H_eq1(k,k))^2/(sum_inf+K/(down_SNR*beta1^2)));
        end
    end
    C1(i_e) = temp/N_iter;
    C2(i_e) = temp1/N_iter;
end
figure('color',[1,1,1,])
plot(NMSE,C1,'k-o',...
   NMSE,C2,'r-^','Linewidth',1.5); 
grid on 
xlabel('NMSE (dB)');
ylabel('Achievable sum-rate (bits/s/Hz)');
L1=legend('Beam selection with perfect CSI','Beam selection with imperfect CSI','Location','NorthEast');set(gca,'FontSize',11, 'FontName','Arial')
set(gca,'XDir','reverse')
% set(gca,'YTick',[20:10:90])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

axis([-30,0,30,100])