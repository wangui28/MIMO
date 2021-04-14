function [h_est, MSE_list] = CRAF_CE(y, W, sys, options, h)
%CRAF0 Summary of this function goes here
%   Detailed explanation goes here
D1 = exp(-1i*2*pi*(0:1:(sys.N1-1))'*(0:1:(options.n1_dic-1))/options.n1_dic)/sqrt(sys.N1);
D2 = exp(-1i*2*pi*(0:1:(sys.N2-1))'*(0:1:(options.n2_dic-1))/options.n2_dic)/sqrt(sys.N2);
D = kron(D1, D2);
Amatrix = W'*D;

options_CRAF.m = options.M;
options_CRAF.n1 = sys.N;
options_CRAF.T = 100;
options_CRAF.cplx_flag = 1;
options_CRAF.betaRAF = 0.1;
options_CRAF.muRAF = 0.8 * (1 - options_CRAF.cplx_flag) + 1 * options_CRAF.cplx_flag;
options_CRAF.tol = 1e-5;
    
[h_est_sparse, MSE_list] = CRAF(abs(y).^2, D'*h, options_CRAF, Amatrix, options.k, 1);
h_est = D*h_est_sparse;

end

