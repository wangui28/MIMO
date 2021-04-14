function [ options ] = MyDefaultOptions(sys)
%MYDEFAULTOPTIONS Summary of this function goes here
% Detailed explanation goes here
options = struct;

options.M = sys.M;  %����ֵ����M=BNrf
options.N_block = sys.M/sys.block_size;      %Ԥ��ʱ��֡��B
% define an over complete dictionary for initialization
if (strcmp(sys.array_type,'UPA'))
    options.n1_dic = sys.N1;
    options.n2_dic = sys.N2;
end
    
options.tol = 1e-8;          %Ԥ����ֹ��ֵ
options.init_supp = 20;         %Ԥ������������
options.k = 20;            % Ԥ��ϡ���on-grid method: how many non-zero elements to estimate?
options.Np = 10;             % Ԥ��·����off-grid method: how many paths to estimate?

end

