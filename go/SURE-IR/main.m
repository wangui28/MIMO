close all;clear all;clc;
n=64;
m=32;
K=8;
PSNR=25;
%---------------------------------------amplitude£º  ap
d=random('unif',0,2*pi,K,1);
ap=exp(1i*d);
%---------------------------------------frequency£º  f0
f0=random('unif',0,1,K,1);
dist=pdist(f0,'euclid');
while min(dist)<1/n
    f0=random('unif',0,1,K,1);
    dist=pdist(f0,'euclid');
end
%---------------------------------------sampling time£º smp
idx=randperm(n);
smp_pmt=idx(1:m);
smp=sort(smp_pmt,'ascend');
%---------------------------------------original signal y_true
y_true=exp(-2i*pi*(1:n)'*f0')*ap;
y=exp(-2i*pi*smp'*f0')*ap;
% --------------------------------------noise measurements y
sigma = (10^(PSNR/10))^-1;
noise_t = sqrt(sigma)/sqrt(2)*randn(n,1) + 1i*sqrt(sigma)/sqrt(2)*randn(n,1);
y = y + noise_t(smp,:);
%---------------------------------------Proposed method
[fes,xes,Rsnr,succ]=SURE_IR(y,n,m,f0,smp,y_true)




