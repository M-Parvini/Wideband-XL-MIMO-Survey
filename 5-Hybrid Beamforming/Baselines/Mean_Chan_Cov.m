function [FRF, FBB, WRF, WBB] = MCCM(H, Chan, OFDM, BS, UE, MO)
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
%%% FRF 
Rt = zeros(Nt, Nt);
for n=1:N
    Rt = Rt + (1/N)*H(:,:,n)'*H(:,:,n);
end
[~,~,Vt] = svd(Rt);
F_MCCM = Vt(:,1:NRF);
FRF = 1/sqrt(Nt)*exp(1i*angle(F_MCCM));

%%% WRF 
Rr = zeros(Nr, Nr);
for n=1:N
    Rr = Rr + (1/N)*H(:,:,n)*H(:,:,n)';
end
[Ur,~,~] = svd(Rr);
W_MCCM = Ur(:,1:NRF);
WRF = 1/sqrt(Nr)*exp(1i*angle(W_MCCM));

%% FBB and WBB
for m = 1:N
    Heq = WRF'*H(:,:,m)*FRF;
    [S,~,V] = svd(Heq);
    FBB(:,:,m) = V(:,1:Ns);
    WBB(:,:,m) = S(:,1:Ns);
end

for m = 1:N
    FBB(:,:,m) = sqrt(Ns) * FBB(:,:,m) / norm(FRF * FBB(:,:,m),'fro');
    WBB(:,:,m) = sqrt(Ns) * WBB(:,:,m) / norm(WRF * WBB(:,:,m),'fro');
end
