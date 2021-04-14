
function [fp_dm_up,power_up]= PSO_Algorithm_up(CUE,DUE,GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C)

%% ��ʼ����Ⱥ
N = 150;                         % ��ʼ��Ⱥ����
d = 2;                          % �ռ�ά��
ger = 120;                      % ����������  
%��ţ��
step0=0.9;                      %��Ⱥ��Χ
step1=0.2;
step=10;%��ʼ����
eta=0.95;
c=2;
k=0.4;
d0=step/c;
limit =zeros(d,2);
limit(1)=0;
limit(2)=0.25 ;              % ����λ�ò�������
vlimit(1)= -0.1;               % �����ٶ�����
vlimit(2)= 0.1 ;
w = 0.8;                        % ����Ȩ��
c1 = 1;                       % ����ѧϰ����
c2 = 1;                       % Ⱥ��ѧϰ����
num=0;
p_c_up=0.25.*rand(N, d);
v = rand(N, d);                  % ��ʼ��Ⱥ���ٶ�
p_cm_up = p_c_up;                          % ÿ���������ʷ���λ��
p_dm_up = zeros(1, d);                % ��Ⱥ����ʷ���λ��
fp_cm_up = -1./zeros(N,1);               % ÿ���������ʷ�����Ӧ��
fp_dm_up = -inf;                      % ��Ⱥ��ʷ�����Ӧ��
%plot(xm,func_pso(xm), 'ro');title('��ʼ״̬ͼ');
%% Ⱥ�����
iter = 1;
record = zeros(ger, 1);          % ��¼��
average=record ;
while iter <= ger
    fp_c_up=zeros(N,1);
    fp_c_up = func_pso_up(p_c_up(:,1),p_c_up(:,2),GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C) ; % ���嵱ǰ��Ӧ��   
     for i = 1:N      
        if fp_cm_up(i) < fp_c_up(i)   %������Ⱥ
            fp_cm_up(i) = fp_c_up(i);     % ���¸�����ʷ�����Ӧ��
            p_cm_up(i,:) = p_c_up(i,:);   % ���¸�����ʷ���λ��
        end 
     end
if fp_dm_up < max(fp_cm_up)    %��ʷ���
        [fp_dm_up , nmax] = max(fp_cm_up);   % ����Ⱥ����ʷ�����Ӧ��
        p_dm_up = p_cm_up(nmax, :);  % ����Ⱥ����ʷ���λ��
 end
    v = v * w + c1 * rand * (p_cm_up - p_c_up) + c2 * rand * (repmat(p_dm_up, N, 1) - p_c_up);% �ٶȸ���
    % �߽��ٶȴ���
    for ii=1:N
        for jj=1:d
            
     if v(ii,jj)>vlimit(2)  v(ii,jj)= vlimit(2);end
      if v(ii,jj)<vlimit(1)  v(ii,jj)= vlimit(1);end
        end
    end
     %BAS����λ���ƶ�
        step=eta*step;
        p_c_left= p_c_up+v.*d0/2;
        fleft=func_pso_up(p_c_left(:,1),p_c_left(:,2),GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C);
        p_c_right=p_c_up-v.*d0/2;
        fright=func_pso_up(p_c_right(:,1),p_c_right(:,2),GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C);
        p=step*v.*sign(fleft-fright);
    p_c_up = p_c_up + k*v+(1-k)*p;% λ�ø���
    % �߽�λ�ô���
        for ii=1:N
        for jj=1:d         
             if p_c_up(ii,jj)>limit(2)  p_c_up(ii,jj)= limit(2);end
              if p_c_up(ii,jj)<limit(1)  p_c_up(ii,jj)= limit(1);end
        end
        end
   
    record(iter) = fp_dm_up;%���ֵ��¼
   average(iter)=mean(fp_cm_up);
%    plot(x, func_pso(x,N), 'ro');title('״̬λ�ñ仯')
     pause(0.1)
    iter = iter+1;
    %disp(['��',num2str(iter-1),'�ε���''���ֵ��',num2str(fym),'����ȡֵ��',num2str(ym)]);
end
% figure(2);
% plot(record,'b');hold on
% plot(average,'r');
% title('��������');
% xlabel('��������');
% ylabel('Ŀ�꺯��ֵ');
% legend('������Ӧ��','ƽ����Ӧ��');
% % figure(3);plot(x, func_pso(x,N), 'ro');title('����״̬λ��');
fp_dm_up;
power_up=sum(p_dm_up);
end
% disp(['��Ч���ֵ��',num2str(fp_dm)]);%����ֵ�����ַ������
% disp(['CUE��DUE�Ĺ���ȡֵ��',num2str(p_dm)]);

