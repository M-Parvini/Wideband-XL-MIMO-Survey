function [OFDM, Chan, BS, UE] = InitializeParams(SNRdBList, Nt, Nr)

OFDM.SNRdB = 0;                 % Signal-to-noise ratio
OFDM.SNRdBList = SNRdBList;     % Signal-to-noise ratio
OFDM.numStreams = 1;            % Number of parallel data streams (Ns)
OFDM.bps = 4;                   % Bits per QAM symbol
OFDM.nfft = 1024;               % FFT length
OFDM.cpLen = 50;                % Cyclic prefix length
OFDM.numOFDMSym = 1;            % Number of OFDM symbols
OFDM.subs = 15e3;               % Subcarrier spacing (B/M)
OFDM.BW = OFDM.nfft*OFDM.subs;  % System Bandwidth
OFDM.NF = 9;                    % Noise figure (not used when SNR is used)

%%% Channel parameters %%%
Chan.delay_spread = 1000e-9;    % Delay Spread of the channel
Chan.doppler = 0;
Chan.fc = 2e9;                  % center frequency
Chan.LSpeed = physconst('LightSpeed');
Chan.lambda = Chan.LSpeed/Chan.fc;
Chan.ChannelType = 'Custom';
Chan.LoS = true;
Chan.LoSKfactor = 1;
Chan.numPaths = 40;              % Will be determined after sampling

Chan.pathDelays = [0 3 7 9 11 19 41]*Chan.delay_spread;   % Path Delays
Chan.pathGains  = [0 -1 -2 -3 -8 -17.2 -20.8];

% [Chan.pathGains, Chan.pathDelays] = Exp_PDP(OFDM.BW, Chan.delay_spread, Chan.numPaths);
%%% Check false for Noise free transmission
Chan.Noise = true;

%%% BS and UE parameters; LoS angle calculation %%%
BS.angle = [0; 0];              % Initial angles (not fixed)
BS.location = [0; 0];           % BS location
BS.AntennaType = 'ULA';         % ULA, UPA
BS.nAntenna = Nt;               % Number of transmit antennas (Nt)
UE.angle = [0; 0];              % Initial angles (not fixed)
UE.location = [1; 1.7321];      % UE location
UE.AntennaType = 'ULA';         % ULA, UPA
UE.nAntenna = Nr;               % Number of receive antennas (Nr)

% LoS Angle: In case we have LoS angle, find the direct path angle between
% the BS and MS
if Chan.LoS
    LoS.angle = atan((UE.location(2)-BS.location(2))/(UE.location(1)-BS.location(1)));
    % setting the BS and UE antenna angles towards LoS angle
    BS.angle = [LoS.angle; 0];
    UE.angle = [LoS.angle+pi; 0];
end
