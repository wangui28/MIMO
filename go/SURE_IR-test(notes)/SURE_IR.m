function [f_es,x_es,RSNR,succ,lambda]=SURE_IR(y,T,M,f0,smp,y_true)
grid_numN = 96; %网格点数
f_grid=0:1/grid_numN:1-1/grid_numN;
f=f_grid; %初始频率网络
threshold=0.005; %修剪阈值
index_amp=1:grid_numN; %索引
gamma = 1; %初始gamma
max_iter=160; %最大迭代次数


for n_itr=1:max_iter
    if norm(z_old-z_new) < 1e-6 %向量2-范数
        break;
    end %判断是否结束迭代，跳出循环
    
    if gamma>1e-8 && norm(z_old-z_new)<gamma^0.5
        gamma=gamma/10; %epsilon十分之一递减
    end %控制epsilon范围，逐步细化
    
    
    z_old=z_new; %更新估计值
    dd=1./(abs(z_old).^2+gamma); 
    D=diag(dd); %生成预测对角阵
    
    %% pruning and lambda update
    lambda = 0.01; %初始正则化参数lambda=100
    if gamma<1e-5
        index_amp = 1:length(index_amp);
        index_t=find(abs(z_new)>threshold); %找出较大系数
        if ~isempty(index_t) %非空矩阵
            index_amp=index_t; %去除较小系数
        end %修剪操作
        
        f=f(index_amp); %去除对应频率分量
        D=D(index_amp,index_amp); %降低对角矩阵维度
        z_old = z_old(index_amp); %去除较小估计量
        A=exp(-2i*pi*smp'*f); %基频率分量矩阵
        lambda = max( 5*(norm(y-A*z_old,'fro')^2 )/M,1e-8); %更新lambda
    end
    
    f_new=f; %更新频率网络
    
    %% sequential update of dictionary parameters
    for i=1:length(index_amp)
        A=exp(-2i*pi*smp'*f); %基频率分量矩阵

        A_d=-2i*pi*exp(-2i*pi*smp'*f);
        dA_theta=zeros(M,length(f));
        dA_theta(:,i)=smp'.*A_d(:,i); %矩阵A对theta求一阶导

        ddtt=-inv(lambda*D+A'*A)*(dA_theta'*A+A'*dA_theta)*inv(lambda*D+A'*A);
        d_obj_theta_s = trace( -(y*y')*( dA_theta*inv(lambda*D+A'*A)*A'+A*ddtt*A'+A*inv(lambda*D+A'*A)*dA_theta' ) );
        der=real(d_obj_theta_s); %f(theta)关于theta的一阶导
        
        stepsize=1;
        f_new(i)=f(i)-stepsize*der;
        f_new=mod(f_new,1); %更新频率
        A_new=exp(-2i*pi*smp'*f_new); %更新基频率分量
        func_val=-y'*A_new*inv(lambda*D+A_new'*A_new)*A_new'*y; %新的代理函数（t+1次）
        sur_val=-y'*A*inv(lambda*D+A'*A)*A'*y; %已知代理函数（t次）
        
        maxit=0;
        while func_val>sur_val && maxit<3
            stepsize=0.1*stepsize; %步长十分之一精度递减
            f_new(i)=f(i)-stepsize*der;
            f_new=mod(f_new,1);
            A_new=exp(-2i*pi*smp'*f_new);
            func_val=-y'*A_new*inv(lambda*D+A_new'*A_new)*A_new'*y;
            maxit=maxit+1;
        end %寻找合适的频率theta满足替代函数递减
        f=f_new; %重建频率网络
    end
    A= exp(-2i*pi*smp'*f);
    z_new=(lambda*D+A'*A)^-1*A'*y;
end

f_es= f; %频率
x_es = z_new; %振幅

succ=0;
if length(f0)==length(f_es)
    if norm(sort(f0)-sort(f_es)') < 1e-3 %确知频率矩阵与重建频率矩阵的误差
        succ=1; 
    end
end %判断算法成功

RSNR= 20*log10(norm(y_true)/norm(y_true-exp(-2i*pi*(1:T)'*f_es)*x_es)); %重构信噪比

end


