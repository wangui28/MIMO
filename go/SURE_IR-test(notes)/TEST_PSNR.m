close all;clear;clc;
T=64; %网格点数 2pi/n
M=32; %测量信号数
K=8; %分辨频率分量数
psnr=[0 5 10 15 20 25 30 35 40]; %峰值信噪比
Sr=zeros(1,9);
Avr=zeros(1,9);
    
for j=2:9
    
    PSNR=psnr(j);
    u=zeros(1,100);
    v=zeros(1,100);

for i=1:100
%---------------------------------------amplitude：  ap
d=random('unif',0,2*pi,K,1); %生成k行1列、[0，2pi]上的均匀分布随机数列
ap=exp(1i*d); %1i虚数单位，生成复振幅
%---------------------------------------frequency：  f0
f0=random('unif',0,1,K,1); %生成k行1列、[0，1]上的均匀分布随机数列
dist=pdist(f0,'euclid'); %矩阵f0中元素的距离
while min(dist)<1/T
    f0=random('unif',0,1,K,1);
    dist=pdist(f0,'euclid');
end %保证频率矩阵中元素距离均大于1/n
%---------------------------------------sample： smp
idx= randperm(T); %将n个数（1-n）随机排序
smp_pmt= idx(1:M); %返回矩阵idx前m个元素
smp= sort(smp_pmt,'ascend'); %将矩阵smp_pmt中元素做升序排列，随机采样
%---------------------------------------original signal： y_true
y_true= exp(-2i*pi*(1:T)'*f0')*ap; %生成确知信号
y_test= exp(-2i*pi*smp'*f0')*ap; %生成抽样测试信号
% --------------------------------------noise measurements: y
sigma= (10^(PSNR/10))^-1; %噪声标准差
noise_t= sqrt(sigma)/sqrt(2)*randn(T,1) + 1i*sqrt(sigma)/sqrt(2)*randn(T,1); %随机产生均值为0的正态分布噪声信号
y= y_test + noise_t(smp,:); %叠加前m个噪声信号
%---------------------------------------Proposed method
[f_es,x_es,RSNR,succ,lambda]= SURE_IR(y,T,M,f0,smp,y_true); %调用算法 [频率 幅度 重构信噪比 成功率 正则化参数]

u(i)=succ;
v(i)=RSNR;
end

sn=sum(u(:)==1);
Sr(j)=sn/100

av=mean(v(:));
Avr(j)=av

end

figure(1);plot(psnr,Sr,'-*b');
set(gca,'XTick',0:5:40); 
set(gca,'YTick',0:0.1:1);
legend('成功率曲线');   
xlabel('峰值信噪比PSNR(dB)');ylabel('成功率');grid on;

figure(2);plot(psnr,Avr,'-*b');
set(gca,'XTick',0:5:40); 
set(gca,'YTick',0:5:60);
legend('重构信噪比曲线');
xlabel('峰值信噪比PSNR(dB)');ylabel('重构信噪比(dB)');grid on;