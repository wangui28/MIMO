 function [theta_es,z_es,err]=IR_SURE_CE(Y,X,W,Nx,Nt,Nr,Ny,Rth)%定义输入和输出
max_outer_iter=3;
max_inner_iter=500;
Snum1 = 5;
Snum2 = 10;
theta = zeros(2,0);%2行的空矩阵
z_old = zeros(1,0);%1行的空矩阵
for outer_iter = 1:max_outer_iter%1，2，3
    At=exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
    Ar=exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
    R = Y-W*Ar*diag(z_old)*At'*X;  %冗余项
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
    theta_update = zeros(2,Snum1);%(2,5)
    for i = 1:Snum1
        u = U(:,i);%取出第i列
        [~,ui] = max(fft(W'*u));
        theta_update(2, i) = (Nr-ui+1)/Nr;%放在第2行
        v = V(:,i);
        [~,vi] = max(fft(X*v));
        theta_update(1, i) = (Nt-vi+1)/Nt;%放在第1行
    end
    theta = [theta theta_update];
    epsilon = 1;
    z_new=[z_old;ones(Snum1,1)];%1行+5行1列
    z_old=[z_old;zeros(Snum1,1)];

    index_amp=1:numel(z_old);

    stepsize_old = 1;
    for inner_itr=1:max_inner_iter%500
        if norm(z_old-z_new)<1e-20
            %break;
            % if |z_old-z_new|<1e-10, break and return to final result
        end
        if epsilon>1e-8 && norm(z_old-z_new)<epsilon^0.5
            epsilon=epsilon/sqrt(10);
        end
        z_old=z_new;
 
        dd=1./(abs(z_old).^2+epsilon);%对应公式（12）
        D=diag(dd);
   
        lambda = 10;
        % pruning and lambda update
        if epsilon<1e-3
            L_index_amp0 = length(index_amp);
            index_amp = 1:L_index_amp0;
            threshold=0.005;
            if (numel(z_new) > Snum2)%如果大于10
                z_sort = sort(abs(z_new),'descend');%降序排列
                threshold = max(z_sort(Snum2+1),0.005);%阈值取第11大的数与0.005中较大值
            end
            index_t=find(abs(z_new)>threshold);%z_new中比阈值大的元素序号
            if ~isempty(index_t)%如果不是空矩阵
                index_amp=index_t;
            end
            % prune the small components of z
            D=D(index_amp,index_amp);%对角线上的元素
            z_old = z_old(index_amp);
            theta = theta(:,index_amp);
            L_index_amp2 = numel(index_amp);
            At=exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
            Ar=exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
            R = Y-W*Ar*diag(z_old)*At'*X;
            Rnorm = norm(R,'fro');
            lambda = max( 1*(Rnorm^2),1e-8);%公式（15）
            % update the weight parameter lambda
        end

        theta_new = theta;
        L_new = length(index_amp);
        %% gradient descend
        dtheta = zeros(2, L_new);
        At = exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
        Ar = exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        At_multiply_X = At'*X;
        At_multiply_X2 = (At_multiply_X*At_multiply_X').';
        Ar_multiply_W = Ar'*W';
        Ar_multiply_W2 = Ar_multiply_W*Ar_multiply_W';
        sigma_ky = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky = sigma_ky + (Ar_multiply_W*Y(:,p)).*conj(At_multiply_X(:,p));
        end
        sigma_kk = At_multiply_X2.*Ar_multiply_W2;
        inv_dkk_multiply_sigma_ky = (D/lambda + sigma_kk)\(sigma_ky);%公式（13）
        pAt = diag(-2i*pi*(0:1:Nt-1))*exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
        pAr = diag(-2i*pi*(0:1:Nr-1))*exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        pAt_multiply_X = pAt'*X;
        pAr_multiply_W = pAr'*W';
        for i = 1:L_new
            p1 = zeros(1,L_new);
            p2 = zeros(1,L_new);
            for ii = 1:Nx
                p1(1,i) = p1(1,i) + Y(:,ii)'*Ar_multiply_W(i,:)'*pAt_multiply_X(i,ii);
                p2(1,i) = p2(1,i) + Y(:,ii)'*pAr_multiply_W(i,:)'*At_multiply_X(i,ii);
            end
            p30 = zeros(L_new,L_new);
            p40 = zeros(L_new,L_new);
            p30(:,i) = At_multiply_X*X'*pAt(:,i);
            p30(i,:) = p30(i,:) + p30(:,i)';
            p40(:,i) = Ar_multiply_W*W*pAr(:,i);
            p40(i,:) = p40(i,:) + p40(:,i)';
            p3 = p30.'.*Ar_multiply_W2;
            p4 = At_multiply_X2.*p40;

            df_dtheta_r_l = -p2*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p2'+inv_dkk_multiply_sigma_ky'*p4*inv_dkk_multiply_sigma_ky;
            df_dtheta_t_l = -p1*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p1'+inv_dkk_multiply_sigma_ky'*p3*inv_dkk_multiply_sigma_ky;

            dtheta(1, i) = real(df_dtheta_t_l);
            dtheta(2, i) = real(df_dtheta_r_l);
        end
        stepsize = stepsize_old*4;%步长变大
        theta_new = theta - stepsize*dtheta;%角度更新
        theta_new = mod(theta_new, 1);

        At_new = exp(-2i*pi*(0:1:Nt-1)'*theta_new(1,:));
        Ar_new = exp(-2i*pi*(0:1:Nr-1)'*theta_new(2,:));
        At_multiply_X_new = At_new'*X;
        At_multiply_X2_new = (At_multiply_X_new*At_multiply_X_new').';
        Ar_multiply_W_new = Ar_new'*W';
        Ar_multiply_W2_new = Ar_multiply_W_new*Ar_multiply_W_new';
        sigma_ky_new = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky_new = sigma_ky_new + (Ar_multiply_W_new*Y(:,p)).*conj(At_multiply_X_new(:,p));
        end
        sigma_kk_new = At_multiply_X2_new.*Ar_multiply_W2_new;

        func_val = -sigma_ky_new'*((D/lambda + sigma_kk_new)\sigma_ky_new);%新
        sur_val = -sigma_ky'*inv_dkk_multiply_sigma_ky;%旧

        maxit=0;
        while func_val>sur_val && maxit<50%新的值大于旧的值
            stepsize = 0.1*stepsize;%减小步长
            theta_new = theta - stepsize*(dtheta/norm(dtheta,'fro'));
            theta_new = mod(theta_new, 1);

            At_new = exp(-2i*pi*(0:1:Nt-1)'*theta_new(1,:));
            Ar_new = exp(-2i*pi*(0:1:Nr-1)'*theta_new(2,:));
            At_multiply_X_new = At_new'*X;
            At_multiply_X2_new = (At_multiply_X_new*At_multiply_X_new').';
            Ar_multiply_W_new = Ar_new'*W';
            Ar_multiply_W2_new = Ar_multiply_W_new*Ar_multiply_W_new';
            sigma_ky_new = zeros(L_new,1);
            for p = 1:Nx
                sigma_ky_new = sigma_ky_new + (Ar_multiply_W_new*Y(:,p)).*conj(At_multiply_X_new(:,p));
            end
            sigma_kk_new = At_multiply_X2_new.*Ar_multiply_W2_new;

            func_val = -sigma_ky_new'*((D/lambda + sigma_kk_new)\sigma_ky_new);%角度更新后的新的目标函数
            maxit=maxit+1;
        end
        stepsize_old = stepsize;
        if (maxit < 30)
            theta = theta_new;
        end

        At = exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
        Ar = exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        At_multiply_X = At'*X;
        At_multiply_X2 = (At_multiply_X*At_multiply_X').';
        Ar_multiply_W = Ar'*W';
        Ar_multiply_W2 = Ar_multiply_W*Ar_multiply_W';
        sigma_ky = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky = sigma_ky + (Ar_multiply_W*Y(:,p)).*conj(At_multiply_X(:,p));
        end
        sigma_kk = At_multiply_X2.*Ar_multiply_W2;
        
        z_new = (D/lambda + sigma_kk)\(sigma_ky);%公式（13）
    end
    z_old = z_new;
end
[z_es z_sort] = sort(z_old,'descend');
theta_es = theta(:,z_sort);
At = exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
Ar = exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
err = norm(Y-W*Ar*diag(z_es)*At'*X,'fro');
end


