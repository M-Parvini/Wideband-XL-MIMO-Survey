function [Rxnoisy, Chan] = ...
          CDLChannelFreqEstimate(OFDM, Chan, H_f, H_fc, wT, wR)

SigPwr = sum(mean(abs(OFDM.txOut).^2));

channelOut = zeros(OFDM.nfft, length(wR));
for n=1:length(wR)
    temp = reshape(squeeze(H_f(n,:,:)), [], OFDM.nfft).*OFDM.txOut.';
    channelOut(:,n) = sum(temp, 1);
end
% Signal reception + Noise addition
if Chan.Noise
    noiseComplexSig = 1/sqrt(2)*(randn(size(channelOut)) + 1i*randn(size(channelOut)));
    noisePwr = SigPwr/(db2pow(OFDM.SNRdB));
    Chan.noisePwr = noisePwr;
else
    noiseComplexSig = 0;
    noisePwr = 0;
end
Rxnoisy = channelOut + noiseComplexSig*sqrt(noisePwr);

% Equal Channel frequency response
heff_f = zeros(OFDM.nfft, 1);
for n=1:OFDM.nfft
    % heff(i) = sqrt(Chan.Nt*Chan.Nr)*wR'*reshape(squeeze(CFR(i,:,:)), Chan.Nt, Chan.Nr).'*wT;
    heff_f(n) = wR'*squeeze(H_f(:,:,n))*wT;
end

heff_fc = zeros(OFDM.nfft, 1);
for n=1:OFDM.nfft
    % heff(i) = sqrt(Chan.Nt*Chan.Nr)*wR'*reshape(squeeze(CFR(i,:,:)), Chan.Nt, Chan.Nr).'*wT;
    heff_fc(n) = wR'*squeeze(H_fc(:,:,n))*wT;
end
Chan.heff = heff_fc;

