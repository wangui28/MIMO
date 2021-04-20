function [ theta1, theta2, alpha ] = NOMP_init(order, N1, N2, W, r, theta1_previous, theta2_previous, stepscale)
%OMP Summary of this function goes here
% Detailed explanation goes here
[~, ~] = size(W);
%P = numel(theta1_previous);
theta1 = theta1_previous;
theta2 = theta2_previous;  %初始值

delta = 1e-8;
G = zeros(3,3);
for i = 1:3
    theta1_i = theta1 + delta*(i-2);
    for j = 1:3
        theta2_j = theta2 + delta*(j-2);
        a = kron(exp(1i*2*pi*(0:1:(N1-1))*theta1_i)', exp(1i*2*pi*(0:1:(N2-1))*theta2_j)');
        f = W'*a;
        fr = f'*r;
        fnorm = f'*f;
        G(i,j) = fr'*fr/fnorm;  
    end
end
I = zeros(2,1);
I(1) = (G(3,1)+G(3,2)+G(3,3)-G(1,1)-G(1,2)-G(1,3))/6/delta;
I(2) = (G(1,3)+G(2,3)+G(3,3)-G(1,1)-G(2,1)-G(3,1))/6/delta;
H = zeros(2,2);
H(1,1) = (G(3,1)+G(3,2)+G(3,3)+G(1,1)+G(1,2)+G(1,3)-2*G(2,1)-2*G(2,2)-2*G(2,3))/3/delta/delta;
H(1,2) = (G(1,1)+G(3,3)-G(1,3)-G(3,1))/4/delta/delta;
H(2,1) = H(1,2);
H(2,2) = (G(1,3)+G(2,3)+G(3,3)+G(1,1)+G(2,1)+G(3,1)-2*G(1,2)-2*G(2,2)-2*G(3,2))/3/delta/delta;

for i = 1:10
%     theta_new = [theta1;theta2]-stepscale*inv(H)*I;
    HH = inv(H);
    ss = 0.5*(HH(1)+HH(2));
    if (strcmp(order,'1st'))
        theta_new = [theta1;theta2]-stepscale*ss*I;
    elseif (strcmp(order,'2nd'))
        theta_new = [theta1;theta2]-stepscale*inv(H)*I;
    end
    a = kron(exp(1i*2*pi*(0:1:(N1-1))*theta_new(1))', exp(1i*2*pi*(0:1:(N2-1))*theta_new(2))');
    f = W'*a;
    fr = f'*r;
    fnorm = f'*f;
    G_new = fr'*fr/fnorm;
    if (G_new >= G(2,2))
        %disp(['stepscale=',num2str(stepscale)]);
        break;
    else
        stepscale = -stepscale/2;  %减小步长
    end
end

theta1 = theta_new(1);
theta2 = theta_new(2);  %新的路径方向


a = kron(exp(1i*2*pi*(0:1:(N1-1))*theta1)', exp(1i*2*pi*(0:1:(N2-1))*theta2)');
f = W'*a;
fr = f'*r;
alpha = fr/norm(f,2)^2;  %新的路径增益

end