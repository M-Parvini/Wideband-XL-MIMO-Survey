function [Ber, eqSym] = ...
    OFDMBasebandDemodulator(Rxnoisy, dataBitsIn, OFDM, Chan, wR, qamTx)

% OFDM demodulation
ofdmDemodIn = Rxnoisy;
CombinedRx = ofdmDemodIn*conj(wR);
Chan.heff = reshape(Chan.heff, [], 1);
% OFDM equalization
eqSym = CombinedRx./Chan.heff;
dataBitsOut = qamdemod(eqSym,2^OFDM.bps,'gray',OutputType='bit', UnitAveragePower=true);
%%% Noise power calculation
% noncombined = qamTx.*Chan.heff;
% noise = noncombined - CombinedRx;
% disp(pow2db(1/mean(abs(noise).^2)))
% BER calculation
% rxBits = int8(dataBitsOut);

%%% total error vector
errVecTotal=abs(dataBitsIn(:)-dataBitsOut(:));
Ber.AllCarriers = nnz(errVecTotal)/length(dataBitsIn(:));

%%% 50% middle carriers error vector
flag = (OFDM.nfft/4)*OFDM.bps;
MidCarriers = flag+1:3*flag;
errVecMiddle=errVecTotal(MidCarriers);
Ber.MiddleCarriers = nnz(errVecMiddle)/length(MidCarriers);

%%% 50% middle carriers error vector
OutCarriers = setdiff(1:OFDM.nfft*OFDM.bps, MidCarriers); 
errVecOut=errVecTotal(OutCarriers);
Ber.OutCarriers = nnz(errVecOut)/length(OutCarriers);
