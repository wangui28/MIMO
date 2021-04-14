function [ Sum_Rate ] = EvaluateSumRate( sys,H,H_est,noise_var)
%EVALUATESUMRATE Summary of this function goes here
%   Detailed explanation goes here
K = size(H,2);

P_A = zeros(sys.N,K);

for k = 1:K
    phase = angle(H_est(:,k));
    P_A(:,k) = 1/sqrt(sys.N)*exp(1i*phase);
end
H_eq = H_est'*P_A;
P_D = H_eq'/(H_eq*H_eq');

N_TS = 64;
P_D = P_D/norm(P_D,'fro')*sqrt(K*sys.N);
P = P_A*P_D;



R = H'*P*sqrt(sys.N/N_TS);
Sum_Rate = 0;
for k = 1:K
    Sum_Rate = Sum_Rate + log2(1+abs(R(k,k))^2/(sum(abs(R(:,k).^2))-abs(R(k,k))^2 + noise_var));
end


end

