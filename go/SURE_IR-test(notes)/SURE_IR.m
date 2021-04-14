function [f_es,x_es,RSNR,succ,lambda]=SURE_IR(y,T,M,f0,smp,y_true)
grid_numN = 96; %�������
f_grid=0:1/grid_numN:1-1/grid_numN;
f=f_grid; %��ʼƵ������
threshold=0.005; %�޼���ֵ
index_amp=1:grid_numN; %����
gamma = 1; %��ʼgamma
max_iter=160; %����������


for n_itr=1:max_iter
    if norm(z_old-z_new) < 1e-6 %����2-����
        break;
    end %�ж��Ƿ��������������ѭ��
    
    if gamma>1e-8 && norm(z_old-z_new)<gamma^0.5
        gamma=gamma/10; %epsilonʮ��֮һ�ݼ�
    end %����epsilon��Χ����ϸ��
    
    
    z_old=z_new; %���¹���ֵ
    dd=1./(abs(z_old).^2+gamma); 
    D=diag(dd); %����Ԥ��Խ���
    
    %% pruning and lambda update
    lambda = 0.01; %��ʼ���򻯲���lambda=100
    if gamma<1e-5
        index_amp = 1:length(index_amp);
        index_t=find(abs(z_new)>threshold); %�ҳ��ϴ�ϵ��
        if ~isempty(index_t) %�ǿվ���
            index_amp=index_t; %ȥ����Сϵ��
        end %�޼�����
        
        f=f(index_amp); %ȥ����ӦƵ�ʷ���
        D=D(index_amp,index_amp); %���ͶԽǾ���ά��
        z_old = z_old(index_amp); %ȥ����С������
        A=exp(-2i*pi*smp'*f); %��Ƶ�ʷ�������
        lambda = max( 5*(norm(y-A*z_old,'fro')^2 )/M,1e-8); %����lambda
    end
    
    f_new=f; %����Ƶ������
    
    %% sequential update of dictionary parameters
    for i=1:length(index_amp)
        A=exp(-2i*pi*smp'*f); %��Ƶ�ʷ�������

        A_d=-2i*pi*exp(-2i*pi*smp'*f);
        dA_theta=zeros(M,length(f));
        dA_theta(:,i)=smp'.*A_d(:,i); %����A��theta��һ�׵�

        ddtt=-inv(lambda*D+A'*A)*(dA_theta'*A+A'*dA_theta)*inv(lambda*D+A'*A);
        d_obj_theta_s = trace( -(y*y')*( dA_theta*inv(lambda*D+A'*A)*A'+A*ddtt*A'+A*inv(lambda*D+A'*A)*dA_theta' ) );
        der=real(d_obj_theta_s); %f(theta)����theta��һ�׵�
        
        stepsize=1;
        f_new(i)=f(i)-stepsize*der;
        f_new=mod(f_new,1); %����Ƶ��
        A_new=exp(-2i*pi*smp'*f_new); %���»�Ƶ�ʷ���
        func_val=-y'*A_new*inv(lambda*D+A_new'*A_new)*A_new'*y; %�µĴ�������t+1�Σ�
        sur_val=-y'*A*inv(lambda*D+A'*A)*A'*y; %��֪��������t�Σ�
        
        maxit=0;
        while func_val>sur_val && maxit<3
            stepsize=0.1*stepsize; %����ʮ��֮һ���ȵݼ�
            f_new(i)=f(i)-stepsize*der;
            f_new=mod(f_new,1);
            A_new=exp(-2i*pi*smp'*f_new);
            func_val=-y'*A_new*inv(lambda*D+A_new'*A_new)*A_new'*y;
            maxit=maxit+1;
        end %Ѱ�Һ��ʵ�Ƶ��theta������������ݼ�
        f=f_new; %�ؽ�Ƶ������
    end
    A= exp(-2i*pi*smp'*f);
    z_new=(lambda*D+A'*A)^-1*A'*y;
end

f_es= f; %Ƶ��
x_es = z_new; %���

succ=0;
if length(f0)==length(f_es)
    if norm(sort(f0)-sort(f_es)') < 1e-3 %ȷ֪Ƶ�ʾ������ؽ�Ƶ�ʾ�������
        succ=1; 
    end
end %�ж��㷨�ɹ�

RSNR= 20*log10(norm(y_true)/norm(y_true-exp(-2i*pi*(1:T)'*f_es)*x_es)); %�ع������

end


