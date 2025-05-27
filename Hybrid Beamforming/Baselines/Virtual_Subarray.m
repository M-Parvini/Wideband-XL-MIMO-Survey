function [FRF, FBB, WRF, WBB] = Virtual_Subarray(H, Chan, OFDM, BS, UE)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Algorithm
% Sorting the channel gains and angles
[~,I] = sort(abs(gains),'descend');
BSAoD = sin(BSAoD(I))/2; % Normalized angles
UEAoA = sin(UEAoA(I))/2; % Normalized angles
gains = gains(I);

%%% BS side
FRF = zeros(Nt, NRF);
% BS Side beam broadening
for l = 1:NRF
    % array response (Eq. 20) g_l
    gl = zeros(Nt, 1);
    % Maximum Beam Squint
    delta_psi = abs((fs/fc)*BSAoD(l));

    % Subarray size
    Ms = SubarrayCal(delta_psi, Nt);
    % Group size
    V = floor(Nt/Ms);
    for v=1:V
        psi_v = BSAoD(l)+ 0.5/Ms + (1/Ms)*(v-1);
        phi_v = (pi/Ms)*(v-1);
        for m=1:Ms
            phaseshift = ((v-1)*Ms+(m-1))*psi_v;
            gl((v-1)*Ms+m, 1) = 1/sqrt(V*Ms)*exp(-1i*2*pi*phaseshift)*...
                                exp(-1i*phi_v);
        end
    end
    FRF(:, l) = gl;
end

% psdmat = FRF'*FRF;
% inv_psdMat = sqrtm(psdmat);
% for n=1:N
%     H_ = conj(H(:,:,n))*FRF/inv_psdMat;
%     [~,~,D] = svd(H_);
%     FBB(:,:,n) = inv_psdMat\D(:,1:Ns);
% end

%%% UE side
WRF = zeros(Nr, NRF);
% UE Side beam broadening
for l = 1:NRF
    % array response (Eq. 20) g_l
    gl = zeros(Nr, 1);
    % Maximum Beam Squint
    delta_psi = abs((fs/fc)*UEAoA(l));

    % Subarray size
    Ms = SubarrayCal(delta_psi, Nr);
    % Group size
    V = floor(Nr/Ms);
    for v=1:V
        psi_v = UEAoA(l)+ 0.5/Ms + (1/Ms)*(v-1);
        phi_v = (pi/Ms)*(v-1);
        for m=1:Ms
            phaseshift = ((v-1)*Ms+(m-1))*psi_v;
            gl((v-1)*Ms+m, 1) = 1/sqrt(V*Ms)*exp(-1i*2*pi*phaseshift)*...
                                exp(-1i*phi_v);
        end
    end
    WRF(:, l) = gl;
end

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


% psdmat = WRF'*WRF;
% for n=1:N
%     H_ = H(:,:,n)'*WRF*inv(sqrtm(psdmat));
%     [~,~,D] = svd(H_);
%     FBB(:,:,n) = D(:,1:Ns);
% end

end

function Subarraysize = SubarrayCal(delta, M)
    %%%% Calculating the quantity of subarraysize
    Ms1 = round((sqrt(1+4*M*delta)-1)/(2*delta));
    Ms2 = round(sqrt(M/delta));
    MsVec = Ms1:Ms2;
    Groupcondition = (1./MsVec).*floor(M./MsVec);
    [~,idx]=find(Groupcondition>=delta);
    if isempty(idx)
        Subarraysize=M;
    else
        Subarraysize = MsVec(max(idx));
    end
end