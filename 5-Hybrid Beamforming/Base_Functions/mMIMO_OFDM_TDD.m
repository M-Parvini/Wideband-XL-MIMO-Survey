function results = mMIMO_OFDM_TDD(OFDM, Chan, BS, UE, SNR)
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

%% Delay-Phase precoding with TTD (Fully connected)
[FRF, FBB, WRF, WBB] = DPP_TTD_Full(H_fc, Chan, OFDM, BS, UE);
DPPTTD_Full.FRF = FRF; DPPTTD_Full.FBB = FBB; DPPTTD_Full.WRF = WRF; DPPTTD_Full.WBB = WBB;

%% Delay-Phase precoding with TTD (Partially connected)
[FRF, FBB, WRF, WBB] = DPP_TTD_Partial(H_fc, Chan, OFDM, BS, UE);
DPPTTD_Partial.FRF = FRF; DPPTTD_Partial.FBB = FBB; DPPTTD_Partial.WRF = WRF; DPPTTD_Partial.WBB = WBB;

% powers
[P_FC, P_AoSA] = PowerConsumption(Chan, OFDM, BS, UE);

%% Spectral Efficiency Calculation
results.TDD_FC = 0;
results.TDD_AoSA = 0;

%% Energy Efficiency Calculation
results.TDD_FC_EE = 0;
results.TDD_AoSA_EE = 0;
%%%
Ns = OFDM.numOFDMSym;
SNR = db2pow(SNR);
for k = 1:OFDM.nfft
    %% Delay Phase Precoding (Fully Connected)
    results.TDD_FC = results.TDD_FC + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(DPPTTD_Full.WRF(:,:,k) * DPPTTD_Full.WBB(:,:,k)) *...
    ChannelMat(:,:,k) * DPPTTD_Full.FRF(:,:,k) * DPPTTD_Full.FBB(:,:,k) * DPPTTD_Full.FBB(:,:,k)' ...
    * DPPTTD_Full.FRF(:,:,k)' * ChannelMat(:,:,k)'*DPPTTD_Full.WRF(:,:,k) * ...
    DPPTTD_Full.WBB(:,:,k)))/OFDM.nfft;
    
    %% Delay Phase Precoding (Partially Connected)
    results.TDD_AoSA = results.TDD_AoSA + ...
    log2(det(eye(Ns) + SNR/Ns * pinv(DPPTTD_Partial.WRF(:,:,k) * DPPTTD_Partial.WBB(:,:,k)) *...
    ChannelMat(:,:,k) * DPPTTD_Partial.FRF(:,:,k) * DPPTTD_Partial.FBB(:,:,k) * DPPTTD_Partial.FBB(:,:,k)' ...
    * DPPTTD_Partial.FRF(:,:,k)' * ChannelMat(:,:,k)'*DPPTTD_Partial.WRF(:,:,k) * ...
    DPPTTD_Partial.WBB(:,:,k)))/OFDM.nfft;

end
%%% Energy efficiency calculation
results.TDD_FC_EE = results.TDD_FC/P_FC;
results.TDD_AoSA_EE = results.TDD_AoSA/P_AoSA;

end