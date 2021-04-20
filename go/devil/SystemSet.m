function [sys] = SystemSet(Nt, Nr, Nx, Ny, L)
%  System parameters
%
sys = struct;

sys.Nt = Nt; % 发射天线数
sys.Nr = Nr; % 接收天线数
sys.N = sys.Nt*sys.Nr;
sys.Nx = Nx; % 发射导频数
sys.Ny = Ny; % 接收导频数
sys.M = sys.Nx*sys.Ny;
sys.L = L; % 路径数

end

