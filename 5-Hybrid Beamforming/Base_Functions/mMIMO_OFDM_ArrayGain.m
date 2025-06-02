function results = mMIMO_OFDM_ArrayGain(OFDM, Chan, BS, UE, SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASSIVE MIMO
% Beam Squint Effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseband Modulator
% Wireless Channel Modeling
[Chan, H_f, H_fc] = mmWaveChannelModeling(Chan, OFDM, BS, UE);

if Chan.Squint
    ChannelMat = H_f;
else
    ChannelMat = H_fc;
end

% Analog Beamforming
% [wT, wR] = Analog_precoder(Chan, BS, UE);

% Digital precoding
[Fopt, Wopt] = Digital_Precoder(ChannelMat, Chan, OFDM, BS, UE);
Dig.Fopt = Fopt; Dig.Wopt = Wopt;

% Manifold optimization 
[FRF, FBB, WRF, WBB] = MO_Precoder(H_fc, Chan, OFDM, BS, UE);
MO.FRF = FRF; MO.FBB = FBB; MO.WRF = WRF; MO.WBB = WBB;

% Beam Squint Aware Manifold optimization 
% [FRF, FBB, WRF, WBB] = BSAMO_Precoder(ChannelMat, Chan, OFDM, BS, UE);
% BSAMO.FRF = FRF; BSAMO.FBB = FBB; BSAMO.WRF = WRF; BSAMO.WBB = WBB;

% Orthogonal Matching Pursuit
[FRF, FBB, WRF, WBB] = OMP_Precoder(H_fc, Chan, OFDM, BS, UE);
OMP.FRF = FRF; OMP.FBB = FBB; OMP.WRF = WRF; OMP.WBB = WBB;

% Beam Squint Aware Orthogonal Matching Pursuit
% [FRF, FBB, WRF, WBB] = BSAOMP_Precoder(ChannelMat, Chan, OFDM, BS, UE);
% BSAOMP.FRF = FRF; BSAOMP.FBB = FBB; BSAOMP.WRF = WRF; BSAOMP.WBB = WBB;

% Delay-Phase precoding with TTD (Fully connected)
% [FRF, FBB, WRF, WBB] = DPP_TTD_Full(H_fc, Chan, OFDM, BS, UE);
% DPPTTD_Full.FRF = FRF; DPPTTD_Full.FBB = FBB; DPPTTD_Full.WRF = WRF; DPPTTD_Full.WBB = WBB;

% Delay-Phase precoding with TTD (Partially connected)
% [FRF, FBB, WRF, WBB] = DPP_TTD_Partial(H_fc, Chan, OFDM, BS, UE);
% DPPTTD_Partial.FRF = FRF; DPPTTD_Partial.FBB = FBB; DPPTTD_Partial.WRF = WRF; DPPTTD_Partial.WBB = WBB;

% Virtual Subarray (Beam Broadening)
[FRF, FBB, WRF, WBB] = Virtual_Subarray(H_fc, Chan, OFDM, BS, UE);
VirArray.FRF = FRF; VirArray.FBB = FBB; VirArray.WRF = WRF; VirArray.WBB = WBB;

%% MCCM
% [FRF, FBB, WRF, WBB] = Mean_Chan_Cov(H_f, Chan, OFDM, BS, UE);
% MCCM.FRF = FRF; MCCM.FBB = FBB; MCCM.WRF = WRF; MCCM.WBB = WBB;

% Spectral Efficiency Calculation
MOArrayGain = zeros(1,OFDM.nfft);
BSAMOArrayGain = zeros(1,OFDM.nfft);
OMPArrayGain = zeros(1,OFDM.nfft);
BSAOMPArrayGain = zeros(1,OFDM.nfft);
DPPTTDArrayGainFull = zeros(1,OFDM.nfft);
DPPTTDArrayGainPartial = zeros(1,OFDM.nfft);
DigitalArrayGain = zeros(1,OFDM.nfft);
AnalogArrayGain = zeros(1,OFDM.nfft);
VirtualArrayGain = zeros(1,OFDM.nfft);
MCCMArrayGain = zeros(1,OFDM.nfft);

for k = 1:OFDM.nfft
    % MO AltMin ArrayGain
    MOArrayGain(1,k) = ...
    MO.WBB(:,:,k)' * MO.WRF' * ChannelMat(:,:,k) * ...
    MO.FRF * MO.FBB(:,:,k);
    
    % BSA MO AltMin ArrayGain
    % BSAMOArrayGain(1,k) = ...
    % BSAMO.WBB(:,:,k)' * BSAMO.WRF' * ChannelMat(:,:,k) * ...
    % BSAMO.FRF * BSAMO.FBB(:,:,k);
    
    % OMP AltMin ArrayGain
    OMPArrayGain(1,k) = ...
    OMP.WBB(:,:,k)' * OMP.WRF' * ChannelMat(:,:,k) * ...
    OMP.FRF * OMP.FBB(:,:,k);

    % BSA OMP AltMin ArrayGain
    % BSAOMPArrayGain(1,k) = ...
    % BSAOMP.WBB(:,:,k)' * BSAOMP.WRF(:,:,k)' * ChannelMat(:,:,k) * ...
    % BSAOMP.FRF(:,:,k) * BSAOMP.FBB(:,:,k);

    % Fully Digital ArrayGain
    DigitalArrayGain(1,k) = ...
    Dig.Wopt(:,:,k)' * ChannelMat(:,:,k)*Dig.Fopt(:,:,k);

    % Analog ArrayGain
    % AnalogArrayGain(1,k) = (wR' * squeeze(ChannelMat(:,:,k)) * wT);
    
    % DPP-TTD ArrayGain (Fully connected)
    % DPPTTDArrayGainFull(1,k) = DPPTTD_Full.FBB(:,:,k)'*DPPTTD_Full.WRF(:,:,k)' * ...
    % ChannelMat(:,:,k) * DPPTTD_Full.FRF(:,:,k)*DPPTTD_Full.FBB(:,:,k);

    % DPP-TTD ArrayGain (Partially connected)
    % DPPTTDArrayGainPartial(1,k) = DPPTTD_Partial.FBB(:,:,k)'*DPPTTD_Partial.WRF(:,:,k)' * ...
    % ChannelMat(:,:,k) * DPPTTD_Partial.FRF(:,:,k)*DPPTTD_Partial.FBB(:,:,k);

    % Virtual Subarray ArrayGain
    VirtualArrayGain(1,k) = ...
    VirArray.WBB(:,:,k)'*VirArray.WRF' * ChannelMat(:,:,k) * ...
    VirArray.FRF*VirArray.WBB(:,:,k);

    % MCCM ArrayGain
    % MCCMArrayGain(1,k) = ...
    % MCCM.WBB(:,:,k)'*MCCM.WRF' * ChannelMat(:,:,k) * ...
    % MCCM.FRF*MCCM.WBB(:,:,k);
end

Nt = BS.nAntenna;
Nr = UE.nAntenna;
%%% Analog
results.AnalogArrayGain = abs(AnalogArrayGain)/(sqrt(Nt*Nr));

%%% Digital
results.DigitalArrayGain = abs(DigitalArrayGain)/(sqrt(Nt*Nr));

%%% DPP-TTD (Fully connected)
results.DPPTTDArrayGainFull = abs(DPPTTDArrayGainFull)/(sqrt(Nt*Nr));

%%% DPP-TTD (Partially connected)
results.DPPTTDArrayGainPartial = abs(DPPTTDArrayGainPartial)/(sqrt(Nt*Nr));

%%% Virtual subarray
results.VirtualArrayGain = abs(VirtualArrayGain)/(sqrt(Nt*Nr));

%%% MCCM
results.MCCMArrayGain = abs(MCCMArrayGain)/(sqrt(Nt*Nr));

%%% OMP
results.OMPArrayGain = abs(OMPArrayGain)/(sqrt(Nt*Nr));
results.BSAOMPArrayGain = abs(BSAOMPArrayGain)/(sqrt(Nt*Nr));

%%% MO
results.MOArrayGain = abs(MOArrayGain)/(sqrt(Nt*Nr));
results.BSAMOArrayGain = abs(BSAMOArrayGain)/(sqrt(Nt*Nr));
end