function [phase_es,h_es,succ,lambda]=IR(y,h0,phase0,phase,W,system)

h_new = ones(system.N,1); %初始 t+1次 估计矩阵
h_old = h0; %初始 t次 估计矩阵
threshold = 0.005;  %修剪阈值
index_amp = 1:system.N;  %索引
gamma = 1; %初始gamma
max_iter = 50; %最大迭代次数


for n_itr = 1:max_iter
    
    disp(['itr=', num2str(n_itr)]);  %窗口显示仿真进度
    
    if norm(h_old-h_new)<1e-6 %向量2-范数
        break;
    end %判断是否结束迭代，跳出循环
    
    if gamma>1e-8 && norm(h_old-h_new)<gamma^0.5
        gamma = gamma/10; %epsilon十分之一递减
    end %控制epsilon范围，逐步细化
    
    
    h_old = h_new; %更新估计值
    dd = 1./(abs(h_old).^2+gamma); 
    D = diag(dd); %生成预测对角阵
    
    %% pruning and lambda update
    lambda = 0.01; %初始正则化参数lambda=100
    if gamma<1e-5
        index_amp = 1:length(index_amp);
        index_t=find(abs(h_old)>threshold); %找出较大系数
        if ~isempty(index_t) %非空矩阵
            index_amp=index_t; %去除较小系数
        end %修剪操作
       
        D=D(index_amp,index_amp); %降低对角矩阵维度
        h_old = h_old(index_amp); %去除较小估计量
        C = W'*exp(2i*pi*diag(phase));
        lambda = max( 5*(norm(y-C*h_old,'fro')^2 )/M,1e-8); %更新lambda
    end
    
    phase_new = phase;
    %% sequential update of dictionary parameters
    for i=1:length(index_amp)

        C = W'*exp(2i*pi*diag(phase));
        
        C_d = W'*2i*pi*exp(2i*pi*phase');
        dC_phase = zeros(system.M,length(phase));
        dC_phase(:,i)=C_d; %矩阵A对theta求一阶导
        
        ddtt=-inv(lambda*D+C'*C)*(dC_phase'*C+C'*dC_phase)*inv(lambda*D+C'*C);
        d_obj_theta_s = trace( -(y*y')*( dC_phase*inv(lambda*D+C'*C)*C'+C*ddtt*C'+C*inv(lambda*D+C'*C)*dC_phase' ) );
        der=real(d_obj_theta_s); %f(theta)关于theta的一阶导
        
        stepsize=1;
        phase_new(i)=phase(i)-stepsize*der;
        phase_new=mod(phase_new,1); %更新频率
        C_new=W'*exp(2i*pi*diag(phase_new)); %更新基频率分量
        func_val=-y'*C_new*inv(lambda*D+C_new'*C_new)*C_new'*y; %新的代理函数（t+1次）
        sur_val=-y'*C*inv(lambda*D+C'*C)*C'*y; %已知代理函数（t次）
        
        maxit=0;
        while func_val>sur_val && maxit<3
            stepsize=0.1*stepsize; %步长十分之一精度递减
            phase_new(i)=phase(i)-stepsize*der;
            phase_new=mod(phase_new,1);
            C_new=W'*exp(2i*pi*diag(phase_new));
            func_val=-y'*C_new*inv(lambda*D+C_new'*C_new)*C_new'*y;
            maxit=maxit+1;
        end %寻找合适的频率theta满足替代函数递减
        phase=phase_new; %重建频率网络
    end
    C = W'*exp(2i*pi*diag(phase));
    h_new=(lambda*D+C'*C)^-1*C'*y;
end

phase_es= phase; %频率
h_es = h_new; %振幅

succ=0;
if length(phase0)==length(phase_es)
    if norm(sort(phase0)-sort(phase_es)') < 1e-3 %确知频率矩阵与重建频率矩阵的误差
        succ=1; 
    end
end %判断算法成功

%RSNR= 20*log10(norm(y_true)/norm(y_true-exp(-2i*pi*(1:T)'*phase_es)*h_es)); %重构信噪比

end


