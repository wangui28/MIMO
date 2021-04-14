function [sys] = SystemSet()
%SYSTEMSET
%   ϵͳ����
sys = struct; % �½��ṹ��
sys.Nt = 32; % ����������
sys.Nr = 16; % ����������
sys.N = sys.Nt*sys.Nr;
sys.Nx = 8; % ���䵼Ƶ��
sys.Ny = 8; % ���յ�Ƶ��
sys.M = sys.Nx*sys.Ny;
sys.L = 3; % ·����
sys.T = 20;
sys.T1 = 5; % number of iterations for NOMP
sys.T2 = 100; % number of iterations for off-grid PC-CPR
sys.tol = 1e-8; %Ԥ����ֹ��ֵ
sys.init_supp = 20; % ��ʼ��h0
sys.k = 10;
sys.block = 16;
sys.size = 4;
end

