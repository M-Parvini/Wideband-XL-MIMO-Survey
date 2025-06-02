function [FRF, FBB, WRF, WBB] = MO_Precoder(H, Chan, OFDM, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% Number of RF chains
NRF = OFDM.RFchain;
% Number of OFDM symbols
Ns = OFDM.numOFDMSym;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% max time delay
DelaySpread = Chan.delay_spread;
% OFDM subcarrier numbers
K = OFDM.nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:K
    if(rank(H(:,:,k))>=Ns)        
        [U,S,V] = svd(H(:,:,k));
        Fopt(:,:,k) = V([1:Nt],[1:Ns]);
        Wopt(:,:,k) = U([1:Nr],[1:Ns]);
    end
end

% Manifold optimization algorithm
[FRF, FBB, ~] = MOAltMin_LowComplexity(Fopt, NRF, 0);
% parfor k = 1:K
for k = 1:K
    FBB(:,:,k) = sqrt(Ns) * FBB(:,:,k) / norm(FRF * FBB(:,:,k),'fro');
end

[WRF, WBB, ~] = MOAltMin_LowComplexity(Wopt, NRF, 0);

end
