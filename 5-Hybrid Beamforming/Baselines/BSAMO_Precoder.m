function [FRF_, FBB_, WRF_, WBB_] = BSAMO_Precoder(H, Chan, OFDM, BS, UE)
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
N = OFDM.nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:N
    if(rank(H(:,:,k))>=Ns)        
        [U,S,V] = svd(H(:,:,k));
        Fopt(:,:,k) = V([1:Nt],[1:Ns]);
        Wopt(:,:,k) = U([1:Nr],[1:Ns]);
    end
end


% Manifold optimization algorithm
for k = 1:N
    [FRF_(:,:,k), FBB(:,:,k), ~] = MOAltMin_LowComplexity(Fopt(:,:,k), NRF, 0);
end

for k = 1:N
    FBB_(:,:,k) = sqrt(Ns) * FBB(:,:,k) / norm(FRF_(:,:,k) * FBB(:,:,k),'fro');
end

for k = 1:N
    [WRF_(:,:,k), WBB_(:,:,k), ~] = MOAltMin_LowComplexity(Wopt(:,:,k), NRF, 0);
end

end
