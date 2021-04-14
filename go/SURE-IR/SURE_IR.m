function [f_es,x_es,Rsnr,succ]=SURE_IR(y,n,m,f0,smp,y_true)
grid_num = 80;
f_grid=0:1/grid_num:1-1/grid_num;
x_new=ones(grid_num,1);
x_old=zeros(grid_num,1);
f=f_grid;

threshold=0.05;
index_amp=1:grid_num;
epsilon = 1;
max_iter=160;

for n_itr=1:max_iter
    if norm(x_old-x_new)<1e-6
        break;
    end
    if epsilon>1e-8 && norm(x_old-x_new)<epsilon^0.5
        epsilon=epsilon/10;
    end
    x_old=x_new;
    
    dd=1./(abs(x_old).^2+epsilon);
    D=diag(dd);
    
    lambda = 0.01;
    % pruning and lambda update
    if epsilon<1e-5
        index_amp = 1:length(index_amp);
        index_t=find(abs(x_new)>threshold);
        if ~isempty(index_t)
            index_amp=index_t;
        end
        
        f=f(index_amp);
        D=D(index_amp,index_amp);
        x_old = x_old(index_amp);
        A=exp(-2i*pi*smp'*f);
        lambda = max( 5*(norm(y-A*x_old,'fro')^2 )/m,1e-8);
    end
    
    f_new=f;
    
    %% sequential update of dictionary parameters
    for i=1:length(index_amp)
        A=exp(-2i*pi*smp'*f);   

        A_d=-2i*pi*exp(-2i*pi*smp'*f);
        dA_theta=zeros(m,length(f));
        dA_theta(:,i)=smp'.*A_d(:,i);

        ddtt=-inv(lambda*D+A'*A)*(dA_theta'*A+A'*dA_theta)*inv(lambda*D+A'*A);
        d_obj_theta_s = trace( -(y*y')*( dA_theta*inv(lambda*D+A'*A)*A'+A*ddtt*A'+A*inv(lambda*D+A'*A)*dA_theta' ) );
        der=real(d_obj_theta_s);
        
        stepsize=1;
        f_new(i)=f(i)-stepsize*der;
        f_new=mod(f_new,1);
        A_new=exp(-2i*pi*smp'*f_new);
        func_val=-y'*A_new*inv(lambda*D+A_new'*A_new)*A_new'*y;
        sur_val=-y'*A*inv(lambda*D+A'*A)*A'*y;
        
        maxit=0;
        while func_val>sur_val && maxit<3
            stepsize=0.1*stepsize;
            f_new(i)=f(i)-stepsize*der;
            f_new=mod(f_new,1);
            A_new=exp(-2i*pi*smp'*f_new);
            func_val=-y'*A_new*inv(lambda*D+A_new'*A_new)*A_new'*y;
            maxit=maxit+1;
        end
        f=f_new;
    end
    A=exp(-2i*pi*smp'*f);
    x_new=(lambda*D+A'*A)^-1*A'*y;
end

f_es=f;
x_es=x_new;

succ=0;
if length(f0)==length(f_es)
    if norm(sort(f0)-sort(f_es)') < 1e-3
        succ=1;
    end
end

Rsnr=20*log10(norm(y_true)/norm(y_true-exp(-2i*pi*(1:n)'*f_es)*x_es));

end



