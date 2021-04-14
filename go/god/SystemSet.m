function [system] = SystemSet()
%SYSTEMSET
%   ϵͳ����

    system = struct;                   % �½��ṹ��
    system.M = 64;                    % ������
    system.N1 = 16;                    % UPA ��
    system.N2 = 16;                    % UPA ��
    system.N = system.N1*system.N2;    % ��������
    system.N_RF = 16;                  % RF����/��λ����
    system.L = 5;                      % ·����
    system.SNR = 10;                   % �����
    system.N_block = system.M/system.N_RF;   %ʱ��֡��B
    system.kappa = 10;                 % LoS/NLoS���ʱ�
                                       % if it is 1, no LoS path
end

