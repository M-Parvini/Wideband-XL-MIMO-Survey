function [Chan, H_f, H_fc] = mmWaveChannelModeling(Chan, OFDM, BS, UE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% OFDM subcarrier numbers
N = OFDM.nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aligning the path gains
% pdpDelay=Chan.pathDelays;
% pdpPowDb=Chan.pathGains;
% 
% % align channel taps to sampling grid
% tapDelayVec=[0:ceil(1e-5+pdpDelay(end)*fs)]/fs;
% pdpPowLin=10.^(pdpPowDb/10);
% activeTap=0;
% idxAligned=zeros(1,length(tapDelayVec)-1);
% 
% for tapIndex=1:length(tapDelayVec)-1
%   ind=find((pdpDelay>=tapDelayVec(tapIndex) & pdpDelay<tapDelayVec(tapIndex+1))==1);
%   if ~isempty(ind)
%     activeTap=activeTap+1;
%     pdpDelayAligned(activeTap)=tapDelayVec(tapIndex);
%     pdpPowAligned(activeTap)=sum(pdpPowLin(ind));
%     idxAligned(:,tapIndex)=tapIndex;
%   end
% end
% Normalizing output channel gains
% pdpPowAligned = reNormalize(pdpPowAligned);
tmax = Chan.delay_spread;
% Chan.numPaths = length(pdpPowAligned);
% L = Chan.numPaths;
L = Chan.numClusters;
%%% Channel Modeling
H_f=zeros(Nr,Nt,N);
H_fc=zeros(Nr,Nt,N);
%%% Random Directions for the channel
Alpha_BS = rand(L,2)*pi/6+pi/6;
Alpha_UE = rand(L,2)*pi/6+pi/6;
%%% Predefined directions for checking the correctness of the simulation
% Alpha_BS = ones(L,1)*pi/3;
% Alpha_UE = ones(L,1)*pi/3+pi;
%%% Random gains (Complex-Guassian)
gains = randn(L,1)+randn(L,1)*1j;
delays = rand(L,1)*tmax;
% gains = 1/sqrt(2)*(randn(L,1)+randn(L,1)*1j);
% delays = pdpDelayAligned;
%%% Predefined gains for checking the correctness of the simulation
% gains = (randn(L,1)+randn(L,1)*1j);
% delays = zeros(1, L);
%%% Channel Frequency Response creation with subcarrier impact
% frequency dependent array response matrix
Atf  = zeros(Nt, L, N);
Atfc = zeros(Nt, L);
Arf  = zeros(Nr, L, N);
Arfc = zeros(Nr, L);

if Chan.LoS
    for k=1:N
        fk = fs/(N)*(k-1-(N-1)/2);
        f=fc+fk;
        %%% LoS path
        % array responses
        Arf(:, 1, k) = ArrayResponse(UE.angle,UE,Chan,Chan.LSpeed/f);
        Atf(:, 1, k) = ArrayResponse(BS.angle,BS,Chan,Chan.LSpeed/f);
        % channel model
        H_f(:,:,k)=H_f(:,:,k)+sqrt(Chan.LoSKfactor)*gains(1)*...
            exp(-1j*2*pi*delays(1)*f)*Arf(:, 1, k)*Atf(:, 1, k)';
        %%% NLoS paths
        for j=2:L
            % array responses
            Arf(:, j, k) = ArrayResponse(Alpha_UE(j,:),UE,Chan,Chan.LSpeed/f);
            Atf(:, j, k) = ArrayResponse(Alpha_BS(j,:),BS,Chan,Chan.LSpeed/f);
            % channel model
            H_f(:,:,k)=H_f(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arf(:, j, k)*Atf(:, j, k)';
        end        
    end
else
    for k=1:N
        fk = fs/(N)*(k-1-(N-1)/2);
        f=fc+fk;
        %%% NLoS paths
        for j=1:L
            % array responses
            Arf(:, j, k) = ArrayResponse(Alpha_UE(j,:),UE,Chan,Chan.LSpeed/f);
            Atf(:, j, k) = ArrayResponse(Alpha_BS(j,:),BS,Chan,Chan.LSpeed/f);
            % Channel model
            H_f(:,:,k)=H_f(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arf(:, j, k)*Atf(:, j, k)';
        end        
    end
end

%%% Channel Frequency Response creation for Center frequency (fc)
if Chan.LoS
    for k=1:N
        fk = fs/(N)*(k-1-(N-1)/2);
        f=fc+fk;
        %%% LoS path
        % array responses
        Arfc(:, 1, k) = ArrayResponse(UE.angle,UE,Chan,Chan.LSpeed/fc);
        Atfc(:, 1, k) = ArrayResponse(BS.angle,BS,Chan,Chan.LSpeed/fc);
        % channel model
        H_fc(:,:,k)=H_fc(:,:,k)+sqrt(Chan.LoSKfactor)*gains(1)*...
            exp(-1j*2*pi*delays(1)*f)*Arfc(:, 1, k)*Atfc(:, 1, k)';
        %%% NLoS paths
        for j=2:L
            % array responses
            Arfc(:, j, k) = ArrayResponse(Alpha_UE(j,:),UE,Chan,Chan.LSpeed/fc);
            Atfc(:, j, k) = ArrayResponse(Alpha_BS(j,:),BS,Chan,Chan.LSpeed/fc);
            % Channel model
            H_fc(:,:,k)=H_fc(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arfc(:, j, k)*Atfc(:, j, k)';
        end        
    end
else
    for k=1:N
        fk = fs/(N)*(k-1-(N-1)/2);
        f=fc+fk;
        %%% NLoS paths
        for j=1:L
            % array responses
            Arfc(:, j, k) = ArrayResponse(Alpha_UE(j,:),UE,Chan,Chan.LSpeed/fc);
            Atfc(:, j, k) = ArrayResponse(Alpha_BS(j,:),BS,Chan,Chan.LSpeed/fc);
            % Channel model
            H_fc(:,:,k)=H_fc(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arfc(:, j, k)*Atfc(:, j, k)';
        end        
    end
end

H_f = H_f.*sqrt(Nt*Nr/L);
H_fc = H_fc.*sqrt(Nt*Nr/L);

%%% Saving the array response
Chan.Atf  = Atf;
Chan.Atfc = Atfc;
Chan.Arf  = Arf;
Chan.Arfc = Arfc;
Chan.normalizedPathGains = gains;

%%% Saving the BS and UE angles
Chan.BSangle = Alpha_BS;
Chan.UEangle = Alpha_UE;

%%% To observe the angular-delay plots
% H_fc = H_f;
%%% Delay-Angular plot
% HH = dftmtx(128)'*squeeze(H_f)*conj(dftmtx(1024));
% figure
% if 0
%     [X,Y] = meshgrid(1:128,1:1024);
%     surf(X.',Y.',abs(HH))
% end

end