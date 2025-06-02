function [OFDM, Chan, BS, UE] = InitializeParams_Rate(SNRdB, SCS, RF, Nt, Nr)

OFDM.SNRdB = SNRdB;             % Signal-to-noise ratio
OFDM.subs = SCS;
OFDM.RFchain = RF;              % Number of RF Chains (NRF)
OFDM.nfft = 64;                 % FFT length
OFDM.numOFDMSym = 3;            % Number of OFDM symbols (Ns)
OFDM.BW = OFDM.nfft*OFDM.subs;  % System Bandwidth

%%% Channel parameters %%%
Chan.delay_spread = 1000e-9;    % Delay Spread of the channel
Chan.doppler = 0;
Chan.fc = 8e9;                  % Center frequency (FR3)
Chan.LSpeed = physconst('LightSpeed');
Chan.lambda = Chan.LSpeed/Chan.fc;
Chan.ChannelType = 'Custom';
Chan.LoS = true;
Chan.LoSKfactor = 1;
Chan.numClusters = 8;
Chan.Rays = 5;
Chan.TTD = false;                 % Apply True Time Delay network concept
Chan.BSTTDsize = 4;               % True Time Delay size
Chan.BSSubarraySize = 4;          % True Time Delay size

Chan.UETTDsize = 4;               % True Time Delay size
Chan.UESubarraySize = 4;          % True Time Delay size
Chan.Squint = true;               % Account for beam Squint

%%% BS and UE parameters; LoS angle calculation %%%
BS.angle = [0; 0];                % Initial angles (not fixed)
BS.location = [0; 0];             % BS location
BS.AntennaType = 'ULA';           % ULA, UPA
BS.nAntenna = Nt;                 % Number of transmit antennas (Nt)
UE.angle = [0; 0];                % Initial angles (not fixed)
UE.location = [1; 1.7321];        % UE location
UE.AntennaType = 'ULA';           % ULA, UPA
UE.nAntenna = Nr;                 % Number of receive antennas (Nr)

% LoS Angle: In case we have LoS angle, find the direct path angle between
% the BS and MS
if Chan.LoS
    LoS.angle = atan((UE.location(2)-BS.location(2))/(UE.location(1)-BS.location(1)));
    % setting the BS and UE antenna angles towards LoS angle
    BS.angle = [LoS.angle; 0];
    UE.angle = [LoS.angle+pi; 0];
end
