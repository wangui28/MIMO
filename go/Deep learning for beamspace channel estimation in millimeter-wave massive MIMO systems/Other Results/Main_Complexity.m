clc;
clear all;
N=32:48:512;
M=N/2;
S=24;
T1=10;
T2=8;
Nc=4;

no_OMP=2*S*M.*N+(2/3*S^3+3/2*S^2+5/6*S)*M+(S^4+2*S^3+S^2)/4
no_AMP=2*T1*M.*N+6*T1*M
no_LAMP=2*T2*M.*N+3*T2*M+3*T2*N
no_GM_LAMP=2*T2*M.*N+3*T2*M+9*Nc*T2*N+7*T2*N
figure('color',[1,1,1])


plot(N,no_OMP,'ms-',...
    N,no_AMP,'g*-',...
    N,no_LAMP,'bd-',...  
    N,no_GM_LAMP,'r^-','LineWidth',1.5);

L1=legend('OMP [11]','AMP [15]','LAMP [29]','Proposed GM-LAMP','Location','NorthWest');set(gca,'FontSize',11, 'FontName','Arial')
set(L1,'FontSize',11, 'FontName','Arial')

xlabel('The number of antennas{\it N} at the BS','FontSize',11,'FontName','Arial')
ylabel('The number of complex multiplications','FontSize',11,'FontName','Arial')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'XTick',[32:48:512])
% set(gca,'YTick',[-23:4:0])
axis([32,512,0,8967056])
