function [FRF_, FBB_, WRF_, WBB_] = BSAOMP_Precoder(H, Chan, OFDM, BS, UE)
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
if Chan.Squint
    At = Chan.Atf;
    Ar = Chan.Arf;
else
    At = Chan.Atfc;
    Ar = Chan.Arfc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:N
    if(rank(H(:,:,k))>=Ns)        
        [U,S,V] = svd(H(:,:,k));
        Fopt(:,:,k) = V([1:Nt],[1:Ns]);
        Wopt(:,:,k) = U([1:Nr],[1:Ns]);
    end
end
for k = 1:N
    [FRF_(:,:,k), FBB_(:,:,k)] = OMPAltMin(Fopt(:,:,k), NRF, At(:,:,k));
end

for k = 1:N
    FBB_(:,:,k) = sqrt(Ns) * FBB_(:,:,k) / norm(FRF_(:,:,k) * FBB_(:,:,k),'fro');
end

for k = 1:N
    [WRF_(:,:,k), WBB_(:,:,k)] = OMPAltMin(Wopt(:,:,k), NRF, Ar(:,:,k));
end

end