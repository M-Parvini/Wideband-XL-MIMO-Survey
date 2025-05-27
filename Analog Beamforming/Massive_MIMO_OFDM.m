function results = Massive_MIMO_OFDM(OFDM, Chan, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASSIVE MIMO
% Beam Squint Effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseband Modulator
% OFDM Modulator
[qamTx, dataBitsIn] = OFDMBasebandModulator(OFDM);

% Tx and Rx beamformers modeling
[wT, wR] = CDLChannelSteeingVector(Chan, BS, UE);

% precoding the transmitted signal
OFDM.txOut = qamTx * wT.';

% Wireless Channel Modeling
[Chan, H_f, H_fc] = mmWaveChannelModeling(Chan, OFDM, BS, UE);
% Simulation loop over the SNRs
results.BerToT = 0;
results.BerMid = 0;
results.BerOut = 0;
for SNRId = 1:length(OFDM.SNRdBList)
    % operating SNR value
    OFDM.SNRdB = OFDM.SNRdBList(SNRId);

    % Signal transmission through the channel + Channel estimation
    [Rxnoisy, Chan] = ...
            CDLChannelFreqEstimate(OFDM, Chan, H_f, H_fc, wT, wR);
    
    % OFDM Demodulator and Equalizer
    [Ber, eqSym] = OFDMBasebandDemodulator(Rxnoisy, dataBitsIn, OFDM, Chan, wR, qamTx);
    
    % Saving the current BER value
    results.BerToT(SNRId) = Ber.AllCarriers;
    results.BerMid(SNRId) = Ber.MiddleCarriers;
    results.BerOut(SNRId) = Ber.OutCarriers;
    
    % Scatter plot of received symbols
    if 0
    scatterplot(eqSym(:))
    end
end
end