function [theta_es,z_es,err]=IR_SURE_CE_UPA(Y,X,W,Nx,Nt0,Nr0,Ny,Rth)
Nt = Nt0(1)*Nt0(2);
Nr = Nr0(1)*Nr0(2);
max_outer_iter=5;
max_inner_iter=300;
Snum1 = 5;
Snum2 = 10;
theta = zeros(4,0);
z_old = zeros(1,0);
for outer_iter = 1:max_outer_iter
    [~,L] = size(theta);
    At = zeros(Nt,L);
    Ar = zeros(Nr,L);
    for l = 1:L
        At(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta(2,l)));
        Ar(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta(4,l)));
    end
    R = Y-W*Ar*diag(z_old)*At'*X;
    Rnorm = norm(R,'fro');
    if ((outer_iter>1) && (Rnorm < Rth))
        break;
        % check if the residue is small enough
        % if so, break and return to the final result
        % else, maybe some paths are missed, such missing paths cannot be found by gradient descend
        %       then do SVD on the residue and find the missing paths and
        %       do IR-based channel estimation again
    end
    [U,~,V] = svd(R);
    % SVD based preconditioning
    % find the missed theta's in the previous iteration
    % add to the theta-list of the current iteration
    theta_update = zeros(4,Snum1);
    for i = 1:Snum1 %´ý²âÊÔdebug
        u = U(:,i);
        wu_mat = reshape(W'*u,Nr0');
        wu_mat_2dfft = fft(fft(wu_mat,[],1),[],2);
        [wu_ma1,wu_mi1] = max(wu_mat_2dfft);
        [~,wu_mi2] = max(wu_ma1);
        theta_update(4, i) = (Nr0(1)-wu_mi1(wu_mi2)+1)/Nr0(1);
        theta_update(3, i) = (Nr0(2)-wu_mi2+1)/Nr0(2);
        v = V(:,i);
        xv_mat = reshape(X*v,Nt0');
        xv_mat_2dfft = fft(fft(xv_mat,[],1),[],2);
        [xv_ma1,xv_mi1] = max(xv_mat_2dfft);
        [~,xv_mi2] = max(xv_ma1);
        theta_update(2, i) = (Nt0(1)-xv_mi1(xv_mi2)+1)/Nt0(1);
        theta_update(1, i) = (Nt0(2)-xv_mi2+1)/Nt0(2);
    end
    theta = [theta theta_update];
    epsilon = 1;
    z_new=[z_old;ones(Snum1,1)];
    z_old=[z_old;zeros(Snum1,1)];

    index_amp=1:numel(z_old);

    stepsize_old = 1;
    for inner_itr=1:max_inner_iter
        if norm(z_old-z_new)<1e-10
            break;
            % if |z_old-z_new|<1e-10, break and return to final result
        end
        if epsilon>1e-8 && norm(z_old-z_new)<epsilon^0.5
            epsilon=epsilon/sqrt(10);
        end
        z_old=z_new;
 
        dd=1./(abs(z_old).^2+epsilon);
        D=diag(dd);
   
        lambda = 10;
        % pruning and lambda update
        if epsilon<1e-3
            L_index_amp0 = length(index_amp);
            index_amp = 1:L_index_amp0;
            threshold=0.005;
            if (numel(z_new) > Snum2)
                z_sort = sort(abs(z_new),'descend');
                threshold = max(z_sort(Snum2+1),0.005);
            end
            index_t=find(abs(z_new)>threshold);
            if ~isempty(index_t)
                index_amp=index_t;
            end
            % prune the small components of z
            D=D(index_amp,index_amp);
            z_old = z_old(index_amp);
            theta = theta(:,index_amp);
            L_index_amp2 = numel(index_amp);
            [~,L] = size(theta);
            At = zeros(Nt,L);
            Ar = zeros(Nr,L);
            for l = 1:L
                At(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta(2,l)));
                Ar(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta(4,l)));
            end
            R = Y-W*Ar*diag(z_old)*At'*X;
            Rnorm = norm(R,'fro');
            lambda = max( 1*(Rnorm^2),1e-8);
            % update the weight parameter lambda
        end

        theta_new = theta;
        L_new = length(index_amp);
        %% gradient descend
        dtheta = zeros(4, L_new);
        [~,L] = size(theta);
        At = zeros(Nt,L);
        Ar = zeros(Nr,L);
        for l = 1:L
            At(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta(2,l)));
            Ar(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta(4,l)));
        end
        At_multiply_X = At'*X;
        At_multiply_X2 = (At_multiply_X*At_multiply_X').';
        Ar_multiply_W = Ar'*W';
        Ar_multiply_W2 = Ar_multiply_W*Ar_multiply_W';
        sigma_ky = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky = sigma_ky + (Ar_multiply_W*Y(:,p)).*conj(At_multiply_X(:,p));
        end
        sigma_kk = At_multiply_X2.*Ar_multiply_W2;
        inv_dkk_multiply_sigma_ky = (D/lambda + sigma_kk)\(sigma_ky);
%         pAt = diag(-2i*pi*(0:1:Nt-1))*exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
%         pAr = diag(-2i*pi*(0:1:Nr-1))*exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        pAt1 = zeros(Nt,L);
        pAt2 = zeros(Nt,L);
        pAr1 = zeros(Nr,L);
        pAr2 = zeros(Nr,L);
        for l = 1:L
            pAt1(:,l) = kron(diag(-2i*pi*(0:1:Nt0(1)-1))*exp(-1i*2*pi*[0:Nt0(1)-1]'*theta(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta(2,l)));
            pAt2(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta(1,l)), diag(-2i*pi*(0:1:Nt0(2)-1))*exp(-1i*2*pi*[0:Nt0(2)-1]'*theta(2,l)));
            pAr1(:,l) = kron(diag(-2i*pi*(0:1:Nr0(1)-1))*exp(-1i*2*pi*[0:Nr0(1)-1]'*theta(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta(4,l)));
            pAr2(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta(3,l)), diag(-2i*pi*(0:1:Nr0(2)-1))*exp(-1i*2*pi*[0:Nr0(2)-1]'*theta(4,l)));
        end
        pAt1_multiply_X = pAt1'*X;
        pAt2_multiply_X = pAt2'*X;
        pAr1_multiply_W = pAr1'*W';
        pAr2_multiply_W = pAr2'*W';
        for i = 1:L_new
            p1 = zeros(1,L_new);
            p2 = zeros(1,L_new);
            p3 = zeros(1,L_new);
            p4 = zeros(1,L_new);
            for ii = 1:Nx
                p1(1,i) = p1(1,i) + Y(:,ii)'*Ar_multiply_W(i,:)'*pAt1_multiply_X(i,ii);
                p2(1,i) = p2(1,i) + Y(:,ii)'*Ar_multiply_W(i,:)'*pAt2_multiply_X(i,ii);
                p3(1,i) = p3(1,i) + Y(:,ii)'*pAr1_multiply_W(i,:)'*At_multiply_X(i,ii);
                p4(1,i) = p4(1,i) + Y(:,ii)'*pAr2_multiply_W(i,:)'*At_multiply_X(i,ii);
            end
            p50 = zeros(L_new,L_new);
            p60 = zeros(L_new,L_new);
            p70 = zeros(L_new,L_new);
            p80 = zeros(L_new,L_new);
            p50(:,i) = At_multiply_X*X'*pAt1(:,i);
            p50(i,:) = p50(i,:) + p50(:,i)';
            p60(:,i) = At_multiply_X*X'*pAt2(:,i);
            p60(i,:) = p60(i,:) + p60(:,i)';
            p70(:,i) = Ar_multiply_W*W*pAr1(:,i);
            p70(i,:) = p70(i,:) + p70(:,i)';
            p80(:,i) = Ar_multiply_W*W*pAr2(:,i);
            p80(i,:) = p80(i,:) + p80(:,i)';
            p5 = p50.'.*Ar_multiply_W2;
            p6 = p60.'.*Ar_multiply_W2;
            p7 = At_multiply_X2.*p70;
            p8 = At_multiply_X2.*p80;

            df_dtheta_t1_l = -p1*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p1'+inv_dkk_multiply_sigma_ky'*p5*inv_dkk_multiply_sigma_ky;
            df_dtheta_t2_l = -p2*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p2'+inv_dkk_multiply_sigma_ky'*p6*inv_dkk_multiply_sigma_ky;
            df_dtheta_r1_l = -p3*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p3'+inv_dkk_multiply_sigma_ky'*p7*inv_dkk_multiply_sigma_ky;
            df_dtheta_r2_l = -p4*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p4'+inv_dkk_multiply_sigma_ky'*p8*inv_dkk_multiply_sigma_ky;
            
            dtheta(1, i) = real(df_dtheta_t1_l);
            dtheta(2, i) = real(df_dtheta_t2_l);
            dtheta(3, i) = real(df_dtheta_r1_l);
            dtheta(4, i) = real(df_dtheta_r2_l);
        end
        stepsize = stepsize_old*4;
        theta_new = theta - stepsize*dtheta;
        theta_new = mod(theta_new, 1);
        
        At_new = zeros(Nt,L_new);
        Ar_new = zeros(Nr,L_new);
        for l = 1:L
            At_new(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta_new(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta_new(2,l)));
            Ar_new(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta_new(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta_new(4,l)));
        end
        At_multiply_X_new = At_new'*X;
        At_multiply_X2_new = (At_multiply_X_new*At_multiply_X_new').';
        Ar_multiply_W_new = Ar_new'*W';
        Ar_multiply_W2_new = Ar_multiply_W_new*Ar_multiply_W_new';
        sigma_ky_new = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky_new = sigma_ky_new + (Ar_multiply_W_new*Y(:,p)).*conj(At_multiply_X_new(:,p));
        end
        sigma_kk_new = At_multiply_X2_new.*Ar_multiply_W2_new;

        func_val = -sigma_ky_new'*((D/lambda + sigma_kk_new)\sigma_ky_new);
        sur_val = -sigma_ky'*inv_dkk_multiply_sigma_ky;

        maxit=0;
        while func_val>sur_val && maxit<50
            stepsize = 0.1*stepsize;
            theta_new = theta - stepsize*dtheta;
            theta_new = mod(theta_new, 1);

            At_new = zeros(Nt,L_new);
            Ar_new = zeros(Nr,L_new);
            for l = 1:L
                At_new(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta_new(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta_new(2,l)));
                Ar_new(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta_new(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta_new(4,l)));
            end
            At_multiply_X_new = At_new'*X;
            At_multiply_X2_new = (At_multiply_X_new*At_multiply_X_new').';
            Ar_multiply_W_new = Ar_new'*W';
            Ar_multiply_W2_new = Ar_multiply_W_new*Ar_multiply_W_new';
            sigma_ky_new = zeros(L_new,1);
            for p = 1:Nx
                sigma_ky_new = sigma_ky_new + (Ar_multiply_W_new*Y(:,p)).*conj(At_multiply_X_new(:,p));
            end
            sigma_kk_new = At_multiply_X2_new.*Ar_multiply_W2_new;

            func_val = -sigma_ky_new'*((D/lambda + sigma_kk_new)\sigma_ky_new);
            maxit=maxit+1;
        end
        stepsize_old = stepsize;
        if (maxit < 30)
            theta = theta_new;
        end

        [~,L] = size(theta);
        At = zeros(Nt,L);
        Ar = zeros(Nr,L);
        for l = 1:L
            At(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta(2,l)));
            Ar(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta(4,l)));
        end
        At_multiply_X = At'*X;
        At_multiply_X2 = (At_multiply_X*At_multiply_X').';
        Ar_multiply_W = Ar'*W';
        Ar_multiply_W2 = Ar_multiply_W*Ar_multiply_W';
        sigma_ky = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky = sigma_ky + (Ar_multiply_W*Y(:,p)).*conj(At_multiply_X(:,p));
        end
        sigma_kk = At_multiply_X2.*Ar_multiply_W2;
        
        z_new = (D/lambda + sigma_kk)\(sigma_ky);
    end
    z_old = z_new;
end
[z_es z_sort] = sort(z_old,'descend');
theta_es = theta(:,z_sort);
[~,L] = size(theta_es);
At = zeros(Nt,L);
Ar = zeros(Nr,L);
for l = 1:L
     At(:,l) = kron(exp(-1i*2*pi*[0:Nt0(1)-1]'*theta_es(1,l)), exp(-1i*2*pi*[0:Nt0(2)-1]'*theta_es(2,l)));
     Ar(:,l) = kron(exp(-1i*2*pi*[0:Nr0(1)-1]'*theta_es(3,l)), exp(-1i*2*pi*[0:Nr0(2)-1]'*theta_es(4,l)));
end
err = norm(Y-W*Ar*diag(z_es)*At'*X,'fro');
end


