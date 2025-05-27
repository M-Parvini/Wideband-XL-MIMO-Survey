function results = mMIMO_OFDM_Compare1(OFDM, Chan, BS, UE, SNR)
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

%% Analog Beamforming
% [wT, wR] = Analog_precoder(Chan, BS, UE);

%% Digital precoding
[Fopt, Wopt] = Digital_Precoder(ChannelMat, Chan, OFDM, BS, UE);
Dig.Fopt = Fopt; Dig.Wopt = Wopt;

% %% Manifold optimization 
[FRF, FBB, WRF, WBB] = MO_Precoder(H_fc, Chan, OFDM, BS, UE);
MO.FRF = FRF; MO.FBB = FBB; MO.WRF = WRF; MO.WBB = WBB;

% Beam Squint Aware Manifold optimization 
[FRF, FBB, WRF, WBB] = BSAMO_Precoder(H_f, Chan, OFDM, BS, UE);
BSAMO.FRF = FRF; BSAMO.FBB = FBB; BSAMO.WRF = WRF; BSAMO.WBB = WBB;

% %% Orthogonal Matching Pursuit
[FRF, FBB, WRF, WBB] = OMP_Precoder(H_fc, Chan, OFDM, BS, UE);
OMP.FRF = FRF; OMP.FBB = FBB; OMP.WRF = WRF; OMP.WBB = WBB;

% Beam Squint Aware Orthogonal Matching Pursuit
[FRF, FBB, WRF, WBB] = BSAOMP_Precoder(H_f, Chan, OFDM, BS, UE);
BSAOMP.FRF = FRF; BSAOMP.FBB = FBB; BSAOMP.WRF = WRF; BSAOMP.WBB = WBB;

%% Virtual Subarray (Beam Broadening)
[FRF, FBB, WRF, WBB] = Virtual_Subarray(H_fc, Chan, OFDM, BS, UE);
VirArray.FRF = FRF; VirArray.FBB = FBB; VirArray.WRF = WRF; VirArray.WBB = WBB;

%% MCCM
[FRF, FBB, WRF, WBB] = Mean_Chan_Cov(H_f, Chan, OFDM, BS, UE);
MCCM.FRF = FRF; MCCM.FBB = FBB; MCCM.WRF = WRF; MCCM.WBB = WBB;

%% Spectral Efficiency Calculation
results.MORate = 0;
results.BSAMORate = 0;
results.OMPRate = 0;
results.BSAOMPRate = 0;
results.DigitalRate = 0;
results.VirtualRate=0;
results.MCCM=0;
%%%
Ns = OFDM.numOFDMSym;
SNR = db2pow(SNR);
for k = 1:OFDM.nfft

    %% Manifold Optimization
    results.MORate = results.MORate + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(MO.WRF * MO.WBB(:,:,k)) * ChannelMat(:,:,k) * ...
    MO.FRF * MO.FBB(:,:,k) * MO.FBB(:,:,k)' * MO.FRF' * ChannelMat(:,:,k)' ...
    * MO.WRF * MO.WBB(:,:,k)))/OFDM.nfft;

    % BSA MO AltMin rate
    results.BSAMORate = results.BSAMORate + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(BSAMO.WRF(:,:,k) * BSAMO.WBB(:,:,k)) * ...
    ChannelMat(:,:,k) * BSAMO.FRF(:,:,k) * BSAMO.FBB(:,:,k) * BSAMO.FBB(:,:,k)' * ...
    BSAMO.FRF(:,:,k)' * ChannelMat(:,:,k)' * BSAMO.WRF(:,:,k) * ...
    BSAMO.WBB(:,:,k)))/OFDM.nfft;

    %% Orthogonal Matching Pursuit
    results.OMPRate = results.OMPRate + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(OMP.WRF * OMP.WBB(:,:,k)) * ChannelMat(:,:,k) * ...
    OMP.FRF * OMP.FBB(:,:,k) * OMP.FBB(:,:,k)' * OMP.FRF' * ChannelMat(:,:,k)' ...
    * OMP.WRF * OMP.WBB(:,:,k)))/OFDM.nfft;

    % BSA OMP AltMin rate
    results.BSAOMPRate = results.BSAOMPRate + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(BSAOMP.WRF(:,:,k) * BSAOMP.WBB(:,:,k)) ...
    * ChannelMat(:,:,k) * BSAOMP.FRF(:,:,k) * BSAOMP.FBB(:,:,k) * ...
    BSAOMP.FBB(:,:,k)' * BSAOMP.FRF(:,:,k)' * ChannelMat(:,:,k)' ...
    * BSAOMP.WRF(:,:,k) * BSAOMP.WBB(:,:,k)))/OFDM.nfft;

    %% Fully Digital
    results.DigitalRate = results.DigitalRate + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(Dig.Wopt(:,:,k)) * ChannelMat(:,:,k) * ...
    Dig.Fopt(:,:,k) * Dig.Fopt(:,:,k)' * ChannelMat(:,:,k)' * Dig.Wopt(:,:,k)))/OFDM.nfft;

    %% Virtual Subarray
    results.VirtualRate = results.VirtualRate + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(VirArray.WRF * VirArray.WBB(:,:,k)) * ...
    ChannelMat(:,:,k) * VirArray.FRF * VirArray.FBB(:,:,k) * ...
    VirArray.FBB(:,:,k)' * VirArray.FRF' * ChannelMat(:,:,k)' ...
    * VirArray.WRF * VirArray.WBB(:,:,k)))/OFDM.nfft;

    %% MCCM
    results.MCCM = results.MCCM + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(MCCM.WRF * MCCM.WBB(:,:,k)) * ChannelMat(:,:,k) * ...
    MCCM.FRF * MCCM.FBB(:,:,k) * MCCM.FBB(:,:,k)' * MCCM.FRF' * ChannelMat(:,:,k)' ...
    * MCCM.WRF * MCCM.WBB(:,:,k)))/OFDM.nfft;

end
end