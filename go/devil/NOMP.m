function [ x_est, theta_t, theta_r, alpha] = NOMP(order, Nt, Nr, C, y, L, th, N_refine)
%OMP Summary of this function goes here
% Detailed explanation goes here
[N, ~] = size(C');
[~, N_column] = size(y);
At = exp(-1i*2*pi*(0:1:(Nt-1))'*(0:1:(Nt*2-1))/Nt/2)/sqrt(Nt); %�Ƕ����뿪����ʸ��
Ar = exp(-1i*2*pi*(0:1:(Nr-1))'*(0:1:(Nr*2-1))/Nr/2)/sqrt(Nr); %�Ƕ��򵽴﷽��ʸ��
D = kron(conj(At), Ar); %�Ƕ����ֵ�
A = C*D;

x_est = zeros(N,N_column);
for column_idx = 1:N_column
    yl = y(:,column_idx);
    xl = zeros(N,1);
    theta_t = [];
    theta_r = [];
    alpha = [];
    for l = 1:L
        r = yl - C*xl;  %����ǰ l-1 ��·���Ե�ǰ·����Ӱ�죬����в�
        if (r'*r<th)
            break;
        end
        
        %% ��2N1*2N2����������Ѱ��ʼ����
        if (l<=L)
            mf = A'*r;
            normA = diag(A'*A);
            [~,idx] = max(abs(mf)./normA);  %�ҳ����ֵ��λ��

            thetat = round((idx-1)/Nr/2)/Nt/2; 
            thetar = mod(idx-1,Nr*2)/Nr/2;  %��ʼ�����ĽǶ�ֵ

            for i=1:N_refine
                [ thetat, thetar, alpha0] = NOMP_init(order, Nt, Nr, C, r, thetat, thetar, 1);  %��Ѱ�ŵ������ĳ�ʼ����
            end
%           a = kron(exp(1i*2*pi*(0:1:(N1-1))*theta1)', exp(1i*2*pi*(0:1:(N2-1))*theta2)');
%           f = W'*a;
%           fr = f'*r;
%           a0 = fr/norm(f,2)^2;
            theta_t = [theta_t; thetat];
            theta_r = [theta_r; thetar];
            alpha = [alpha; alpha0];
        end
       
%         P = numel(alpha);
%         for zz = 1:0
%             for p = 1:P
%                 xp = zeros(N, 1);
%                 for pp = 1:P
%                     if (pp~=p)
%                         xp = xp + kron(exp(1i*2*pi*(0:1:(N1-1))*om1(pp))', exp(1i*2*pi*(0:1:(N2-1))*om2(pp))')*alpha(pp);
%                     end
%                 end
%                 r = yl - W'*xp;
%                 [ om1(p), om2(p), alpha(p) ] = NOMP_2D_once(N1, N2, W, r, om1(p), om2(p), 1);
%             end
%         end

        P = numel(alpha);
        xl = zeros(N, 1);
        for p = 1:P
            xl = xl + kron(exp(1i*2*pi*(0:1:(Nt-1))*theta_t(p))', exp(1i*2*pi*(0:1:(Nr-1))*theta_r(p))')*alpha(p);
        end  %��õ�ǰ�ŵ��ķ��������
    end
  
    x_est(:,column_idx) = xl;
end

end

