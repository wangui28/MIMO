function  [xhat,dxdr,dxdr1] = eta_gm(rhat,rvar,theta)
K=4; % the number of Gaussian components
P=exp(theta(1,:))./sum(exp(theta(1,:))); % the probability of Gaussian components
mu=theta(2,:);
var=exp(theta(3,:));
U=0; V=0; dU=0; dV=0; dU1=0; dV1=0;
for k=1:K
    tilde_mu_k_x=(mu(k)*rvar+rhat*var(k))./(rvar+var(k));
    CNpdf_k=exp(-conj(rhat-mu(k)).*(rhat-mu(k))./2./(rvar+var(k)))./2./pi./(var(k)+rvar);
    u_k=tilde_mu_k_x.*P(k).*CNpdf_k;
    du_k=var(k)./(rvar+var(k)).*P(k).*CNpdf_k+u_k.*(-conj(rhat-mu(k))./2./(rvar+var(k)));
    du_k1=u_k.*(-(rhat-mu(k))./2./(rvar+var(k)));
    v_k=P(k).*CNpdf_k;
    dv_k=v_k.*(-conj(rhat-mu(k))./2./(rvar+var(k)));
    dv_k1=v_k.*(-(rhat-mu(k))./2./(rvar+var(k)));
    U = U + u_k;
    dU = dU + du_k;
    dU1 = dU1 +du_k1;
    V = V + v_k;
    dV = dV + dv_k;
    dV1 = dV1 +dv_k1;
end
xhat=U./V;
dxdr=(dU.*V-U.*dV)./(V.*conj(V));
dxdr1=(dU1.*V-U.*dV1)./(V.*conj(V));
dxdr=mean(dxdr,1);
dxdr1=mean(dxdr1,1);
end

