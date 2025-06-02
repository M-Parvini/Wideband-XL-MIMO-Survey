function [FRF, FBB, WRF, WBB] = DPP_TTD_Full(H, Chan, OFDM, BS, UE)
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
% OFDM subcarrier numbers
N = OFDM.nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Required AoA, AoD and Path gains
BSAoD = Chan.BSangle;
UEAoA = Chan.UEangle;
gains = Chan.normalizedPathGains;
%%%%%%%%%%%%%%%%%%%%%%%%% Channel sorting %%%%%%%%%%%%%%%%%%%%%%%%%
[~,I] = sort(abs(gains),'descend');
BSAoD = BSAoD(I); % Normalized angles
UEAoA = UEAoA(I); % Normalized angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BS side TTD correction %%%%%%%%%%%%%%%%%%%%%%%%%
K = Chan.BSTTDsize;
Pt = Nt/K;
BSAoD  = BSAoD(1:NRF);
AAtc = Chan.Atfc(:, :, 1);
AAt = AAtc(:, I);
AAt = AAt(:, 1:NRF);

%%% TTD matrices
Au = zeros(Nt, K*NRF);
AmTTD = zeros(K*NRF, NRF, N);
%%%
for nRF = 1:NRF
    WTx = reshape(AAt(:, nRF), [], K);
    WTx_ = mat2cell(WTx, size(WTx,1), ones(1,size(WTx,2)));
    Au(:, (nRF-1)*K+1:nRF*K) = blkdiag(WTx_{:});
end

for m = 1:N
    f = fs/(N)*(m-1-(N-1)/2);
    fm=fc+f;
    epsim = fm/fc;
    TTDCell = {};
    for nRF = 1:NRF
        beta = (epsim-1)*Pt*sin(BSAoD(nRF));
        TTDCell{nRF} = exp(-1i*pi*[0:K-1]*beta).';
    end
    AmTTD(:,:,m) = blkdiag(TTDCell{:});
end
%%%%%%%%%%%%%%%%%%%%%%% UE side TTD correction %%%%%%%%%%%%%%%%%%%%%%%%%
K = Chan.UETTDsize;
Pr = Nr/K;
UEAoD  = UEAoA(1:NRF);
AArc = Chan.Arfc(:, :, 1);
AAr = AArc(:, I);
AAr = AAr(:, 1:NRF);

%%% TTD matrices
Wu = zeros(Nr, K*NRF);
WmTTD = zeros(K*NRF, NRF, N);
%%%
for nRF = 1:NRF
    WTx = reshape(AAr(:, nRF), [], K);
    WTx_ = mat2cell(WTx, size(WTx,1), ones(1,size(WTx,2)));
    Wu(:, (nRF-1)*K+1:nRF*K) = blkdiag(WTx_{:});
end

for m = 1:N
    f = fs/(N)*(m-1-(N-1)/2);
    fm=fc+f;
    epsim = fm/fc;
    TTDCell = {};
    for nRF = 1:NRF
        beta = (epsim-1)*Pr*sin(UEAoD(nRF));
        TTDCell{nRF} = exp(-1i*pi*[0:K-1]*beta).';
    end
    WmTTD(:,:,m) = blkdiag(TTDCell{:});
end

for m = 1:N
    Heq = WmTTD(:,:,m)'*Wu'*H(:,:,m)*Au*AmTTD(:,:,m);
    [S,~,V] = svd(Heq);
    FBB(:,:,m) = V(:,1:Ns);
    WBB(:,:,m) = S(:,1:Ns);
end

for m = 1:N
    FRF(:,:,m) = Au*AmTTD(:,:,m);
    WRF(:,:,m) = Wu*WmTTD(:,:,m);
    FBB(:,:,m) = sqrt(Ns) * FBB(:,:,m) / norm(FRF(:,:,m) * FBB(:,:,m),'fro');
    WBB(:,:,m) = sqrt(Ns) * WBB(:,:,m) / norm(WRF(:,:,m) * WBB(:,:,m),'fro');
end



% function [FRF, FBB, WRF, WBB] = DPP_TTD_Full(H, Chan, OFDM, BS, UE)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Transmit antenna number N
% Nt = BS.nAntenna;
% % Received antenna number M
% Nr = UE.nAntenna;
% % Number of RF chains
% NRF = OFDM.RFchain;
% % Number of OFDM symbols
% Ns = OFDM.numOFDMSym;
% % central frequency fc
% fc = Chan.fc;
% % bandwidth
% fs = OFDM.BW;
% % OFDM subcarrier numbers
% N = OFDM.nfft;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Required AoA, AoD and Path gains
% BSAoD = Chan.BSangle;
% UEAoA = Chan.UEangle;
% gains = Chan.normalizedPathGains;
% %%%%%%%%%%%%%%%%%%%%%%%%% Channel sorting %%%%%%%%%%%%%%%%%%%%%%%%%
% [~,I] = sort(abs(gains),'descend');
% BSAoD = BSAoD(I); % Normalized angles
% UEAoA = UEAoA(I); % Normalized angles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%% BS side TTD correction %%%%%%%%%%%%%%%%%%%%%%%%%
% K = Chan.BSTTDsize;
% Pt = Nt/K;
% BSAoD  = BSAoD(1:1);
% AAtc = Chan.Atfc(:, :, 1);
% AAt = AAtc(:, I);
% AAt = AAt(:, 1:1);
% 
% %%% TTD matrices
% Au = zeros(Nt, K*NRF);
% AmTTD = zeros(K*NRF, NRF, N);
% %%%
% for nRF = 1:NRF
%     WTx = reshape(AAt(:, 1), [], K);
%     WTx_ = mat2cell(WTx, size(WTx,1), ones(1,size(WTx,2)));
%     Au(:, (nRF-1)*K+1:nRF*K) = blkdiag(WTx_{:});
% end
% 
% for m = 1:N
%     f = fs/(N)*(m-1-(N-1)/2);
%     fm=fc+f;
%     epsim = fm/fc;
%     TTDCell = {};
%     for nRF = 1:NRF
%         beta = (epsim-1)*Pt*sin(BSAoD(1));
%         TTDCell{nRF} = exp(-1i*pi*[0:K-1]*beta).';
%     end
%     AmTTD(:,:,m) = blkdiag(TTDCell{:});
% end
% %%%%%%%%%%%%%%%%%%%%%%% UE side TTD correction %%%%%%%%%%%%%%%%%%%%%%%%%
% K = Chan.UETTDsize;
% Pr = Nr/K;
% UEAoD  = UEAoA(1:1);
% AArc = Chan.Arfc(:, :, 1);
% AAr = AArc(:, I);
% AAr = AAr(:, 1:1);
% 
% %%% TTD matrices
% Wu = zeros(Nr, K*NRF);
% WmTTD = zeros(K*NRF, NRF, N);
% %%%
% for nRF = 1:NRF
%     WTx = reshape(AAr(:, 1), [], K);
%     WTx_ = mat2cell(WTx, size(WTx,1), ones(1,size(WTx,2)));
%     Wu(:, (nRF-1)*K+1:nRF*K) = blkdiag(WTx_{:});
% end
% 
% for m = 1:N
%     f = fs/(N)*(m-1-(N-1)/2);
%     fm=fc+f;
%     epsim = fm/fc;
%     TTDCell = {};
%     for nRF = 1:NRF
%         beta = (epsim-1)*Pr*sin(UEAoD(1));
%         TTDCell{nRF} = exp(-1i*pi*[0:K-1]*beta).';
%     end
%     WmTTD(:,:,m) = blkdiag(TTDCell{:});
% end
% 
% for m = 1:N
%     Heq = WmTTD(:,:,m)'*Wu'*H(:,:,m)*Au*AmTTD(:,:,m);
%     [S,~,V] = svd(Heq);
%     FBB(:,:,m) = V(:,1:Ns);
%     WBB(:,:,m) = S(:,1:Ns);
% end
% 
% for m = 1:N
%     FRF(:,:,m) = Au*AmTTD(:,:,m);
%     WRF(:,:,m) = Wu*WmTTD(:,:,m);
%     FBB(:,:,m) = sqrt(Ns) * FBB(:,:,m) / norm(FRF(:,:,m) * FBB(:,:,m),'fro');
%     WBB(:,:,m) = sqrt(Ns) * WBB(:,:,m) / norm(WRF(:,:,m) * WBB(:,:,m),'fro');
% end
