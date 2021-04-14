function [ x_est, om1, om2, alpha ] = NOMP_2D(order, N1, N2, W, y, L, th, N_refine)
%OMP Summary of this function goes here
% Detailed explanation goes here
[N, ~] = size(W);
[~, N_column] = size(y);
D1 = exp(-1i*2*pi*(0:1:(N1-1))'*(0:1:(N1*2-1))/N1/2)/sqrt(N1);
D2 = exp(-1i*2*pi*(0:1:(N2-1))'*(0:1:(N2*2-1))/N2/2)/sqrt(N2);
D = kron(D1, D2);
WD = W'*D;
x_est = zeros(N,N_column);
for column_idx = 1:N_column
    yl = y(:,column_idx);
    xl = zeros(N,1);
    om1 = [];
    om2 = [];
    alpha = [];
    for l = 1:L
        r = yl - W'*xl;  %消除前 l-1 条路经对当前路径的影响，计算残差
        if (r'*r<th)
            break;
        end
        
        %% 在2N1*2N2的网格中搜寻初始估计
        if (l<=L)
            mf = WD'*r;
            normWD = diag(WD'*WD);
            [~,idx] = max(abs(mf)./normWD);  %找出最大值的位置

            theta1 = round((idx-1)/N2/2)/N1/2; 
            theta2 = mod(idx-1,N2*2)/N2/2;  %初始给定的角度值

            for i=1:N_refine
                [ theta1, theta2, a0 ] = NOMP_2D_once(order, N1, N2, W, r, theta1, theta2, 1);  %搜寻信道参数的初始估计
            end
%           a = kron(exp(1i*2*pi*(0:1:(N1-1))*theta1)', exp(1i*2*pi*(0:1:(N2-1))*theta2)');
%           f = W'*a;
%           fr = f'*r;
%           a0 = fr/norm(f,2)^2;
            om1 = [om1; theta1];
            om2 = [om2; theta2];
            alpha = [alpha; a0];
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
            xl = xl + kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))')*alpha(p);
        end  %获得当前信道的非网格估计
    end
    
    
    x_est(:,column_idx) = xl;
end

end

