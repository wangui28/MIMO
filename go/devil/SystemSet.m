function [sys] = SystemSet()
%SYSTEMSET
%   系统参数
sys = struct; % 新建结构体
sys.Nt = 32; % 发射天线数
sys.Nr = 16; % 接收天线数
sys.N = sys.Nt*sys.Nr;
sys.Nx = 8; % 发射导频数
sys.Ny = 8; % 接收导频数
sys.M = sys.Nx*sys.Ny;
sys.L = 3; % 路径数
sys.T = 20;
sys.T1 = 5; % number of iterations for NOMP
sys.T2 = 100; % number of iterations for off-grid PC-CPR
sys.tol = 1e-8; %预置终止阈值
sys.init_supp = 20; % 初始化h0
sys.k = 10;
sys.block = 16;
sys.size = 4;
end

