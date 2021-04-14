function [system] = SystemSet()
%SYSTEMSET
%   系统参数

    system = struct;                   % 新建结构体
    system.M = 64;                    % 测量数
    system.N1 = 16;                    % UPA 行
    system.N2 = 16;                    % UPA 列
    system.N = system.N1*system.N2;    % 天线总数
    system.N_RF = 16;                  % RF链数/相位块数
    system.L = 5;                      % 路径数
    system.SNR = 10;                   % 信噪比
    system.N_block = system.M/system.N_RF;   %时间帧数B
    system.kappa = 10;                 % LoS/NLoS功率比
                                       % if it is 1, no LoS path
end

