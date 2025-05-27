function [FRF, FBB, WRF, WBB] = DPP_TTD_Partial(H, Chan, OFDM, BS, UE)
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
% BS subarraysize
BS_Subarray_size = Chan.BSSubarraySize;
% UE subarraysize
UE_Subarray_size = Chan.UESubarraySize;
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
Nsub = Nt/NRF;
Ksub = Chan.BSTTDsize;
NTTD = Nsub/Ksub;

BSAoD  = BSAoD(1:NRF);
AAtc = Chan.Atfc(:, :, 1);
AAt = AAtc(:, I);
AAt = AAt(:, 1:NRF);

%%% TTD matrices
FTTDm = zeros(Ksub*NRF, NRF, N);
%%%

temp = cell(1, Ksub*NRF);
for nRF = 1:NRF
    WTx = reshape(AAt((nRF-1)*Nsub+1:nRF*Nsub, nRF), [], Ksub);
    WTx_ = mat2cell(WTx, size(WTx,1), ones(1,size(WTx,2)));
    temp((nRF-1)*Ksub+1:nRF*Ksub) = WTx_;
end

FPS = blkdiag(temp{:});

for m = 1:N
    f = fs/(N)*(m-1-(N-1)/2);
    fm=fc+f;
    epsim = fm/fc;
    TTDCell = {};
    for nRF = 1:NRF
        beta = NTTD*(epsim-1)*sin(BSAoD(nRF));
        TTDCell{nRF} = exp(-1i*pi*(0:Ksub-1)*beta).';
    end
    FTTDm(:,:,m) = blkdiag(TTDCell{:});
end
%%%%%%%%%%%%%%%%%%%%%%% UE side TTD correction %%%%%%%%%%%%%%%%%%%%%%%%%
Nsub = Nr/NRF;
Ksub = Chan.UETTDsize;
NTTD = Nsub/Ksub;

UEAoA  = UEAoA(1:NRF);
AArc = Chan.Arfc(:, :, 1);
AAr = AArc(:, I);
AAr = AAr(:, 1:NRF);

%%% TTD matrices
WTTDm = zeros(Ksub*NRF, NRF, N);
%%%

temp = cell(1, Ksub*NRF);
for nRF = 1:NRF
    WRx = reshape(AAr((nRF-1)*Nsub+1:nRF*Nsub, nRF), [], Ksub);
    WRx_ = mat2cell(WRx, size(WRx,1), ones(1,size(WRx,2)));
    temp((nRF-1)*Ksub+1:nRF*Ksub) = WRx_;
end

WPS = blkdiag(temp{:});

for m = 1:N
    f = fs/(N)*(m-1-(N-1)/2);
    fm=fc+f;
    epsim = fm/fc;
    TTDCell = {};
    for nRF = 1:NRF
        beta = NTTD*(epsim-1)*sin(UEAoA(nRF));
        TTDCell{nRF} = exp(-1i*pi*[0:Ksub-1]*beta).';
    end
    WTTDm(:,:,m) = blkdiag(TTDCell{:});
end

for m = 1:N
    Heq = WTTDm(:,:,m)'*WPS'*H(:,:,m)*FPS*FTTDm(:,:,m);
    [S,~,V] = svd(Heq);
    FBB(:,:,m) = V(:,1:Ns);
    WBB(:,:,m) = S(:,1:Ns);
end

for m = 1:N
    FRF(:,:,m) = FPS*FTTDm(:,:,m);
    WRF(:,:,m) = WPS*WTTDm(:,:,m);
    FBB(:,:,m) = sqrt(Ns) * FBB(:,:,m) / norm(FRF(:,:,m) * FBB(:,:,m),'fro');
    WBB(:,:,m) = sqrt(Ns) * WBB(:,:,m) / norm(WRF(:,:,m) * WBB(:,:,m),'fro');
end


% function [FRF, FBB, WRF, WBB] = DPP_TTD_Partial(H, Chan, OFDM, BS, UE)
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
% % BS subarraysize
% BS_Subarray_size = Chan.BSSubarraySize;
% % UE subarraysize
% UE_Subarray_size = Chan.UESubarraySize;
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
% Nsub = Nt/NRF;
% Ksub = Chan.BSTTDsize;
% NTTD = Nsub/Ksub;
% 
% BSAoD  = BSAoD(1);
% AAtc = Chan.Atfc(:, :, 1);
% AAt = AAtc(:, I);
% AAt = AAt(:, 1);
% 
% %%% TTD matrices
% FTTDm = zeros(Ksub*NRF, NRF, N);
% %%%
% 
% temp = cell(1, Ksub*NRF);
% for nRF = 1:NRF
%     WTx = reshape(AAt((nRF-1)*Nsub+1:nRF*Nsub, 1), [], Ksub);
%     WTx_ = mat2cell(WTx, size(WTx,1), ones(1,size(WTx,2)));
%     temp((nRF-1)*Ksub+1:nRF*Ksub) = WTx_;
% end
% 
% FPS = blkdiag(temp{:});
% 
% for m = 1:N
%     f = fs/(N)*(m-1-(N-1)/2);
%     fm=fc+f;
%     epsim = fm/fc;
%     TTDCell = {};
%     for nRF = 1:NRF
%         beta = NTTD*(epsim-1)*sin(BSAoD(1));
%         TTDCell{nRF} = exp(-1i*pi*(0:Ksub-1)*beta).';
%     end
%     FTTDm(:,:,m) = blkdiag(TTDCell{:});
% end
% %%%%%%%%%%%%%%%%%%%%%%% UE side TTD correction %%%%%%%%%%%%%%%%%%%%%%%%%
% Nsub = Nr/NRF;
% Ksub = Chan.UETTDsize;
% NTTD = Nsub/Ksub;
% 
% UEAoA  = UEAoA(1:1);
% AArc = Chan.Arfc(:, :, 1);
% AAr = AArc(:, I);
% AAr = AAr(:, 1:1);
% 
% %%% TTD matrices
% WTTDm = zeros(Ksub*NRF, NRF, N);
% %%%
% 
% temp = cell(1, Ksub*NRF);
% for nRF = 1:NRF
%     WRx = reshape(AAr((nRF-1)*Nsub+1:nRF*Nsub, 1), [], Ksub);
%     WRx_ = mat2cell(WRx, size(WRx,1), ones(1,size(WRx,2)));
%     temp((nRF-1)*Ksub+1:nRF*Ksub) = WRx_;
% end
% 
% WPS = blkdiag(temp{:});
% 
% for m = 1:N
%     f = fs/(N)*(m-1-(N-1)/2);
%     fm=fc+f;
%     epsim = fm/fc;
%     TTDCell = {};
%     for nRF = 1:NRF
%         beta = NTTD*(epsim-1)*sin(UEAoA(1));
%         TTDCell{nRF} = exp(-1i*pi*[0:Ksub-1]*beta).';
%     end
%     WTTDm(:,:,m) = blkdiag(TTDCell{:});
% end
% 
% for m = 1:N
%     Heq = WTTDm(:,:,m)'*WPS'*H(:,:,m)*FPS*FTTDm(:,:,m);
%     [S,~,V] = svd(Heq);
%     FBB(:,:,m) = V(:,1:Ns);
%     WBB(:,:,m) = S(:,1:Ns);
% end
% 
% for m = 1:N
%     FRF(:,:,m) = FPS*FTTDm(:,:,m);
%     WRF(:,:,m) = WPS*WTTDm(:,:,m);
%     FBB(:,:,m) = sqrt(Ns) * FBB(:,:,m) / norm(FRF(:,:,m) * FBB(:,:,m),'fro');
%     WBB(:,:,m) = sqrt(Ns) * WBB(:,:,m) / norm(WRF(:,:,m) * WBB(:,:,m),'fro');
% end
