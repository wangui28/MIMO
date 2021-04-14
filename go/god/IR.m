function [phase_es,h_es,succ,lambda]=IR(y,h0,phase0,phase,W,system)

h_new = ones(system.N,1); %��ʼ t+1�� ���ƾ���
h_old = h0; %��ʼ t�� ���ƾ���
threshold = 0.005;  %�޼���ֵ
index_amp = 1:system.N;  %����
gamma = 1; %��ʼgamma
max_iter = 50; %����������


for n_itr = 1:max_iter
    
    disp(['itr=', num2str(n_itr)]);  %������ʾ�������
    
    if norm(h_old-h_new)<1e-6 %����2-����
        break;
    end %�ж��Ƿ��������������ѭ��
    
    if gamma>1e-8 && norm(h_old-h_new)<gamma^0.5
        gamma = gamma/10; %epsilonʮ��֮һ�ݼ�
    end %����epsilon��Χ����ϸ��
    
    
    h_old = h_new; %���¹���ֵ
    dd = 1./(abs(h_old).^2+gamma); 
    D = diag(dd); %����Ԥ��Խ���
    
    %% pruning and lambda update
    lambda = 0.01; %��ʼ���򻯲���lambda=100
    if gamma<1e-5
        index_amp = 1:length(index_amp);
        index_t=find(abs(h_old)>threshold); %�ҳ��ϴ�ϵ��
        if ~isempty(index_t) %�ǿվ���
            index_amp=index_t; %ȥ����Сϵ��
        end %�޼�����
       
        D=D(index_amp,index_amp); %���ͶԽǾ���ά��
        h_old = h_old(index_amp); %ȥ����С������
        C = W'*exp(2i*pi*diag(phase));
        lambda = max( 5*(norm(y-C*h_old,'fro')^2 )/M,1e-8); %����lambda
    end
    
    phase_new = phase;
    %% sequential update of dictionary parameters
    for i=1:length(index_amp)

        C = W'*exp(2i*pi*diag(phase));
        
        C_d = W'*2i*pi*exp(2i*pi*phase');
        dC_phase = zeros(system.M,length(phase));
        dC_phase(:,i)=C_d; %����A��theta��һ�׵�
        
        ddtt=-inv(lambda*D+C'*C)*(dC_phase'*C+C'*dC_phase)*inv(lambda*D+C'*C);
        d_obj_theta_s = trace( -(y*y')*( dC_phase*inv(lambda*D+C'*C)*C'+C*ddtt*C'+C*inv(lambda*D+C'*C)*dC_phase' ) );
        der=real(d_obj_theta_s); %f(theta)����theta��һ�׵�
        
        stepsize=1;
        phase_new(i)=phase(i)-stepsize*der;
        phase_new=mod(phase_new,1); %����Ƶ��
        C_new=W'*exp(2i*pi*diag(phase_new)); %���»�Ƶ�ʷ���
        func_val=-y'*C_new*inv(lambda*D+C_new'*C_new)*C_new'*y; %�µĴ�������t+1�Σ�
        sur_val=-y'*C*inv(lambda*D+C'*C)*C'*y; %��֪��������t�Σ�
        
        maxit=0;
        while func_val>sur_val && maxit<3
            stepsize=0.1*stepsize; %����ʮ��֮һ���ȵݼ�
            phase_new(i)=phase(i)-stepsize*der;
            phase_new=mod(phase_new,1);
            C_new=W'*exp(2i*pi*diag(phase_new));
            func_val=-y'*C_new*inv(lambda*D+C_new'*C_new)*C_new'*y;
            maxit=maxit+1;
        end %Ѱ�Һ��ʵ�Ƶ��theta������������ݼ�
        phase=phase_new; %�ؽ�Ƶ������
    end
    C = W'*exp(2i*pi*diag(phase));
    h_new=(lambda*D+C'*C)^-1*C'*y;
end

phase_es= phase; %Ƶ��
h_es = h_new; %���

succ=0;
if length(phase0)==length(phase_es)
    if norm(sort(phase0)-sort(phase_es)') < 1e-3 %ȷ֪Ƶ�ʾ������ؽ�Ƶ�ʾ�������
        succ=1; 
    end
end %�ж��㷨�ɹ�

%RSNR= 20*log10(norm(y_true)/norm(y_true-exp(-2i*pi*(1:T)'*phase_es)*h_es)); %�ع������

end


