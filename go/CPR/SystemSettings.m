function [ sys ] = SystemSettings(M, N_RF, SNR)
%SYSTEMSETTINGS Summary of this function goes here
%   Detailed explanation goes here
    sys = struct;                   % create an empty struct
    sys.array_type = 'UPA';         % the type of antenna array, currently support 'UPA' only
    
    sys.M = M;                      % the number of measurements
    sys.N1 = 32;                    % UPA size in dim 1
    sys.N2 = 16;                    % UPA size in dim 2
    sys.N = sys.N1*sys.N2;          % total number of antennas
    sys.block_size = N_RF;          % number of RF chains /size of phase structure block
    sys.L = 5;                      % number of paths
    sys.SNR = SNR;                  % signal-to-noise ratio (SNR)
    sys.N_trial = 2;            	% number of Monte-Carlo trials
    sys.T = 20;                     % number of iterations for on-grid PC-CPR
    sys.T1 = 5;                     % number of iterations for NOMP
    sys.T2 = 100;                   % number of iterations for off-grid PC-CPR
    sys.kappa = 10;                 % The power ratio for LoS by NLoS
                                    % if it is 1, no LoS path
                                    
end

