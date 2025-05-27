function [qamTx, dataBitsIn] = OFDMBasebandModulator(OFDM)
% Baseband modulator consisting of QAM and OFDM modulation

%%% Create data bits
dataBitsIn = randi([0,1],[OFDM.nfft*OFDM.bps OFDM.numOFDMSym ...
    OFDM.numStreams]);

%%% QAM mode
QAM_M = 2^OFDM.bps; % Modulation order
qamTx = qammod(dataBitsIn,QAM_M,'gray',...
    InputType="bit", ...
    UnitAveragePower=true);

%%% OFDM Modulation
% ofdmOut = ofdmmod(qamTx,OFDM.nfft,OFDM.cpLen)*sqrt(OFDM.nfft);