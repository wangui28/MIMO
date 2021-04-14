
function [fp_dm_up,power_up]= PSO_Algorithm_up(CUE,DUE,GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C)

%% 初始化种群
N = 150;                         % 初始种群个数
d = 2;                          % 空间维数
ger = 120;                      % 最大迭代次数  
%天牛须
step0=0.9;                      %种群范围
step1=0.2;
step=10;%初始步长
eta=0.95;
c=2;
k=0.4;
d0=step/c;
limit =zeros(d,2);
limit(1)=0;
limit(2)=0.25 ;              % 设置位置参数限制
vlimit(1)= -0.1;               % 设置速度限制
vlimit(2)= 0.1 ;
w = 0.8;                        % 惯性权重
c1 = 1;                       % 自我学习因子
c2 = 1;                       % 群体学习因子
num=0;
p_c_up=0.25.*rand(N, d);
v = rand(N, d);                  % 初始种群的速度
p_cm_up = p_c_up;                          % 每个个体的历史最佳位置
p_dm_up = zeros(1, d);                % 种群的历史最佳位置
fp_cm_up = -1./zeros(N,1);               % 每个个体的历史最佳适应度
fp_dm_up = -inf;                      % 种群历史最佳适应度
%plot(xm,func_pso(xm), 'ro');title('初始状态图');
%% 群体更新
iter = 1;
record = zeros(ger, 1);          % 记录器
average=record ;
while iter <= ger
    fp_c_up=zeros(N,1);
    fp_c_up = func_pso_up(p_c_up(:,1),p_c_up(:,2),GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C) ; % 个体当前适应度   
     for i = 1:N      
        if fp_cm_up(i) < fp_c_up(i)   %更新种群
            fp_cm_up(i) = fp_c_up(i);     % 更新个体历史最佳适应度
            p_cm_up(i,:) = p_c_up(i,:);   % 更新个体历史最佳位置
        end 
     end
if fp_dm_up < max(fp_cm_up)    %历史最佳
        [fp_dm_up , nmax] = max(fp_cm_up);   % 更新群体历史最佳适应度
        p_dm_up = p_cm_up(nmax, :);  % 更新群体历史最佳位置
 end
    v = v * w + c1 * rand * (p_cm_up - p_c_up) + c2 * rand * (repmat(p_dm_up, N, 1) - p_c_up);% 速度更新
    % 边界速度处理
    for ii=1:N
        for jj=1:d
            
     if v(ii,jj)>vlimit(2)  v(ii,jj)= vlimit(2);end
      if v(ii,jj)<vlimit(1)  v(ii,jj)= vlimit(1);end
        end
    end
     %BAS部分位置移动
        step=eta*step;
        p_c_left= p_c_up+v.*d0/2;
        fleft=func_pso_up(p_c_left(:,1),p_c_left(:,2),GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C);
        p_c_right=p_c_up-v.*d0/2;
        fright=func_pso_up(p_c_right(:,1),p_c_right(:,2),GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C);
        p=step*v.*sign(fleft-fright);
    p_c_up = p_c_up + k*v+(1-k)*p;% 位置更新
    % 边界位置处理
        for ii=1:N
        for jj=1:d         
             if p_c_up(ii,jj)>limit(2)  p_c_up(ii,jj)= limit(2);end
              if p_c_up(ii,jj)<limit(1)  p_c_up(ii,jj)= limit(1);end
        end
        end
   
    record(iter) = fp_dm_up;%最大值记录
   average(iter)=mean(fp_cm_up);
%    plot(x, func_pso(x,N), 'ro');title('状态位置变化')
     pause(0.1)
    iter = iter+1;
    %disp(['第',num2str(iter-1),'次迭代''最大值：',num2str(fym),'变量取值：',num2str(ym)]);
end
% figure(2);
% plot(record,'b');hold on
% plot(average,'r');
% title('收敛过程');
% xlabel('进化代数');
% ylabel('目标函数值');
% legend('最优适应度','平均适应度');
% % figure(3);plot(x, func_pso(x,N), 'ro');title('最终状态位置');
fp_dm_up;
power_up=sum(p_dm_up);
end
% disp(['能效最大值：',num2str(fp_dm)]);%把数值换成字符串输出
% disp(['CUE和DUE的功率取值：',num2str(p_dm)]);

