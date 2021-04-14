f=0:0.1:1;
Sr=[ 0    0.79   0.82600    0.90200   0.918  0.945  0.954   0.9856    0.994500    1.0000    1.0000   ]; 

figure(1);plot(f,Sr,'-*b');
set(gca,'XTick',0:0.1:1); 
set(gca,'YTick',0:0.1:1);
legend('成功率曲线');   
xlabel('频率间隔系数\Deltaf');ylabel('成功率');grid on;

