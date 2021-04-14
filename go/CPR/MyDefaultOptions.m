function [ options ] = MyDefaultOptions(sys)
%MYDEFAULTOPTIONS Summary of this function goes here
% Detailed explanation goes here
options = struct;

options.M = sys.M;  %测量值总数M=BNrf
options.N_block = sys.M/sys.block_size;      %预置时间帧数B
% define an over complete dictionary for initialization
if (strcmp(sys.array_type,'UPA'))
    options.n1_dic = sys.N1;
    options.n2_dic = sys.N2;
end
    
options.tol = 1e-8;          %预置终止阈值
options.init_supp = 20;         %预置最大迭代次数
options.k = 20;            % 预置稀疏度on-grid method: how many non-zero elements to estimate?
options.Np = 10;             % 预置路径数off-grid method: how many paths to estimate?

end

