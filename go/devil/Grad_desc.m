function [ h_new, om1_new, om2_new, alpha_new ] = Grad_desc(order, N1, N2, W, y, om1, om2, alpha, h, stepsize)
%OMP Summary of this function goes here
%   Detailed explanation goes here
[N, ~] = size(W);
P = numel(om1);
residue = y - W'*h;  %?в?

if (strcmp(order, '1st'))
    dloss_dom1 = zeros(P,1);
    dloss_dom2 = zeros(P,1);
    dloss_dalpha = zeros(P,1);
    for p = 1:P
        ap = kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        Wa = W'*ap;
        dWa_dom1 = W'*kron((1i*2*pi*(0:1:(N1-1)))'.*exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        dWa_dom2 = W'*kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', (1i*2*pi*(0:1:(N2-1)))'.*exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        
        dloss_dalpha(p) = -Wa'*residue;
        dloss_dom1(p) = -2*real(alpha(p)*residue'*dWa_dom1);
        dloss_dom2(p) = -2*real(alpha(p)*residue'*dWa_dom2);
    end

    loss = norm(residue);
    stepscale = 1;
    for i = 1:100
        alpha_new = alpha - stepsize * stepscale * conj(dloss_dalpha);
        om1_new = om1 - stepsize * stepscale * dloss_dom1;
        om2_new = om2 - stepsize * stepscale * dloss_dom2;
        h_new = zeros(N, 1);
        for p = 1:numel(om1)
            h_new = h_new + kron(exp(1i*2*pi*(0:1:(N1-1))*om1_new(p))', exp(1i*2*pi*(0:1:(N2-1))*om2_new(p))')*alpha_new(p);
        end
        loss_new = norm(y - W' * h_new);
        if (loss_new < loss)
            break;
        else
            stepscale = 0.5*stepscale;
        end
    end
    
    
%     alpha_new = alpha - stepsize * conj(dloss_dalpha);
%     om1_new = om1 - stepsize * dloss_dom1;
%     om2_new = om2 - stepsize * dloss_dom2;
% 
%     h_new = zeros(N, 1);
%     for p = 1:numel(om1)
%         h_new = h_new + kron(exp(1i*2*pi*(0:1:(N1-1))*om1_new(p))', exp(1i*2*pi*(0:1:(N2-1))*om2_new(p))')*alpha_new(p);
%     end
elseif (strcmp(order, '2nd'))
    for p = 1:P
        % interference cancelation
        h_p = h - kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))')*alpha(p);
        y_p = y - W'*h_p;
        
        % calculate the 1st order derivatives
        Wa = W'*kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        dloss_dalpha = conj(alpha(p))*Wa'*Wa - y_p'*Wa;

        dWa_dom1 = W'*kron((-1i*2*pi*(0:1:(N1-1)))'.*exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        dloss_dom1 = alpha(p)'*alpha(p)*(dWa_dom1'*Wa + Wa'*dWa_dom1) - 2*real(y_p'*alpha(p)*dWa_dom1);

        dWa_dom2 = W'*kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', (-1i*2*pi*(0:1:(N2-1)))'.*exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        dloss_dom2 = alpha(p)'*alpha(p)*(dWa_dom2'*Wa + Wa'*dWa_dom2) - 2*real(y_p'*alpha(p)*dWa_dom2);
        
        % calculate the 2nd order derivatives
        d2loss_dalpha2 = Wa'*Wa;
        
        d2Wa_dom12 = W'*kron((-1i*2*pi*(0:1:(N1-1)))'.*(-1i*2*pi*(0:1:(N1-1)))'.*exp(1i*2*pi*(0:1:(N1-1))*om1(p))', exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        d2loss_dom12 = alpha(p)'*alpha(p)*2*real(d2Wa_dom12'*Wa + dWa_dom1'*dWa_dom1) - 2*real(y_p'*alpha(p)*d2Wa_dom12);
        
        d2Wa_dom22 = W'*kron(exp(1i*2*pi*(0:1:(N1-1))*om1(p))', (-1i*2*pi*(0:1:(N2-1)))'.*(-1i*2*pi*(0:1:(N2-1)))'.*exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        d2loss_dom22 = alpha(p)'*alpha(p)*2*real(d2Wa_dom22'*Wa + dWa_dom2'*dWa_dom2) - 2*real(y_p'*alpha(p)*d2Wa_dom22);
        
        d2loss_dalphadom1 = conj(alpha(p))*(dWa_dom1'*Wa + Wa'*dWa_dom1) - y_p'*dWa_dom1;
        d2loss_dalphadom2 = conj(alpha(p))*(dWa_dom2'*Wa + Wa'*dWa_dom2) - y_p'*dWa_dom2;
        d2Wa_dom1dom2 = W'*kron((-1i*2*pi*(0:1:(N1-1)))'.*exp(1i*2*pi*(0:1:(N1-1))*om1(p))', (-1i*2*pi*(0:1:(N2-1)))'.*exp(1i*2*pi*(0:1:(N2-1))*om2(p))');
        d2loss_dom1dom2 = alpha(p)'*alpha(p)*2*real(d2Wa_dom1dom2'*Wa + dWa_dom1'*dWa_dom2) - 2*real(y_p'*alpha(p)*d2Wa_dom1dom2);
        
        % build 1st order derivative vector and 2nd order derivative matrix
        I = [dloss_dalpha;dloss_dom1;dloss_dom2];
        H = [d2loss_dalpha2,d2loss_dalphadom1,d2loss_dalphadom2;d2loss_dalphadom1',d2loss_dom12,d2loss_dom1dom2;d2loss_dalphadom2',d2loss_dom1dom2',d2loss_dom22];
        
        loss = norm(y - W'*h);
        stepscale = 0.1;
        for i = 1:10
            refine = [alpha(p);om1(p);om2(p)] - stepscale*conj(inv(H)*I);
            loss_new = norm(y_p - W' * kron(exp(1i*2*pi*(0:1:(N1-1))*refine(2))', exp(1i*2*pi*(0:1:(N2-1))*refine(3))')*refine(1));
            if (loss_new < loss)
                break;
            else
                stepscale = -0.5*stepscale;
            end
        end
        %disp(['stepscale=',num2str(stepscale),'  loss_new=',num2str(loss_new)]);
            
        alpha(p) = refine(1);
        om1(p) = refine(2);
        om2(p) = refine(3);
        
        h = zeros(N, 1);
        for pp = 1:numel(om1)
            h = h + kron(exp(1i*2*pi*(0:1:(N1-1))*om1(pp))', exp(1i*2*pi*(0:1:(N2-1))*om2(pp))')*alpha(pp);
        end
    end
    
    h_new = h;
    om1_new = om1;
    om2_new = om2;
    alpha_new = alpha;
end
    
%disp(norm(y - W'*h)-norm(y - W'*h_new));

end