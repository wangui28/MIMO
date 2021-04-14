close all;clear;clc;
T=64; %������� 2pi/n
M=32; %�����ź���
K=8; %�ֱ�Ƶ�ʷ�����
psnr=[0 5 10 15 20 25 30 35 40]; %��ֵ�����
Sr=zeros(1,9);
Avr=zeros(1,9);
    
for j=2:9
    
    PSNR=psnr(j);
    u=zeros(1,100);
    v=zeros(1,100);

for i=1:100
%---------------------------------------amplitude��  ap
d=random('unif',0,2*pi,K,1); %����k��1�С�[0��2pi]�ϵľ��ȷֲ��������
ap=exp(1i*d); %1i������λ�����ɸ����
%---------------------------------------frequency��  f0
f0=random('unif',0,1,K,1); %����k��1�С�[0��1]�ϵľ��ȷֲ��������
dist=pdist(f0,'euclid'); %����f0��Ԫ�صľ���
while min(dist)<1/T
    f0=random('unif',0,1,K,1);
    dist=pdist(f0,'euclid');
end %��֤Ƶ�ʾ�����Ԫ�ؾ��������1/n
%---------------------------------------sample�� smp
idx= randperm(T); %��n������1-n���������
smp_pmt= idx(1:M); %���ؾ���idxǰm��Ԫ��
smp= sort(smp_pmt,'ascend'); %������smp_pmt��Ԫ�����������У��������
%---------------------------------------original signal�� y_true
y_true= exp(-2i*pi*(1:T)'*f0')*ap; %����ȷ֪�ź�
y_test= exp(-2i*pi*smp'*f0')*ap; %���ɳ��������ź�
% --------------------------------------noise measurements: y
sigma= (10^(PSNR/10))^-1; %������׼��
noise_t= sqrt(sigma)/sqrt(2)*randn(T,1) + 1i*sqrt(sigma)/sqrt(2)*randn(T,1); %���������ֵΪ0����̬�ֲ������ź�
y= y_test + noise_t(smp,:); %����ǰm�������ź�
%---------------------------------------Proposed method
[f_es,x_es,RSNR,succ,lambda]= SURE_IR(y,T,M,f0,smp,y_true); %�����㷨 [Ƶ�� ���� �ع������ �ɹ��� ���򻯲���]

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
legend('�ɹ�������');   
xlabel('��ֵ�����PSNR(dB)');ylabel('�ɹ���');grid on;

figure(2);plot(psnr,Avr,'-*b');
set(gca,'XTick',0:5:40); 
set(gca,'YTick',0:5:60);
legend('�ع����������');
xlabel('��ֵ�����PSNR(dB)');ylabel('�ع������(dB)');grid on;