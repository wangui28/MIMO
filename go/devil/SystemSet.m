function [sys] = SystemSet(Nt, Nr, Nx, Ny, L)
%  System parameters
%
sys = struct;

sys.Nt = Nt; % ����������
sys.Nr = Nr; % ����������
sys.N = sys.Nt*sys.Nr;
sys.Nx = Nx; % ���䵼Ƶ��
sys.Ny = Ny; % ���յ�Ƶ��
sys.M = sys.Nx*sys.Ny;
sys.L = L; % ·����

end

