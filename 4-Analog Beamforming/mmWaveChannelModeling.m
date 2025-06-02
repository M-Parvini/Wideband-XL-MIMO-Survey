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
% max time delay
DelaySpread = Chan.delay_spread;
% OFDM subcarrier numbers
N = OFDM.nfft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aligning the path gains
pdpDelay=Chan.pathDelays;
pdpPowDb=Chan.pathGains;

% align channel taps to sampling grid
tapDelayVec=[0:ceil(1e-5+pdpDelay(end)*fs)]/fs;
pdpPowLin=10.^(pdpPowDb/10);
activeTap=0;
idxAligned=zeros(1,length(tapDelayVec)-1);

for tapIndex=1:length(tapDelayVec)-1
  ind=find((pdpDelay>=tapDelayVec(tapIndex) & pdpDelay<tapDelayVec(tapIndex+1))==1);
  if ~isempty(ind)
    activeTap=activeTap+1;
    pdpDelayAligned(activeTap)=tapDelayVec(tapIndex);
    pdpPowAligned(activeTap)=sum(pdpPowLin(ind));
    idxAligned(:,tapIndex)=tapIndex;
  end
end
% Normalizing output channel gains
pdpPowAligned = reNormalize(pdpPowAligned);
Chan.numPaths = length(pdpPowAligned);
L = Chan.numPaths;
%%% Channel Modeling
H_f=zeros(Nr,Nt,N);
H_fc=zeros(Nr,Nt,N);
%%% Random Directions for the channel
Alpha_BS = rand(L,1)*pi-pi/2;
Alpha_UE = rand(L,1)*pi-pi/2;
%%% Predefined directions for checking the correctness of the simulation
Alpha_BS = ones(L,1)*pi/3;
Alpha_UE = ones(L,1)*pi/3+pi;
%%% Random gains (Complex-Guassian)
% gains = 1/sqrt(2)*(randn(L,1)+randn(L,1)*1j).*sqrt(pdpPowAligned).';
gains = 1/sqrt(2)*(randn(L,1)+randn(L,1)*1j);
delays = pdpDelayAligned;
%%% Predefined gains for checking the correctness of the simulation
% gains = (randn(L,1)+randn(L,1)*1j);
% delays = zeros(1, L);
%%% Channel Frequency Response creation with subcarrier impact
if Chan.LoS
    for k=1:N
        f=fc+fs/(N)*(k-1-(N-1)/2);
        %%% LoS path
        H_f(:,:,k)=H_f(:,:,k)+sqrt(Chan.LoSKfactor)*gains(1)*...
            exp(-1j*2*pi*delays(1)*f)*...
            ArrayResponse(UE.angle,UE,Chan,Chan.LSpeed/f)*...
            ArrayResponse(BS.angle,BS,Chan,Chan.LSpeed/f)';
        %%% NLoS paths
        for j=2:L
            H_f(:,:,k)=H_f(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/f)*...
            ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/f)';
        end        
    end
else
    for k=1:N
        f=fc+fs/(N)*(k-1-(N-1)/2);
        %%% NLoS paths
        for j=1:L
            H_f(:,:,k)=H_f(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/f)*...
            ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/f)';
        end        
    end
end

%%% Channel Frequency Response creation for Center frequency (fc)
if Chan.LoS
    for k=1:N
        f=fc+fs/(N)*(k-1-(N-1)/2);
        %%% LoS path
        H_fc(:,:,k)=H_fc(:,:,k)+sqrt(Chan.LoSKfactor)*gains(1)*...
            exp(-1j*2*pi*delays(1)*f)*...
            ArrayResponse(UE.angle,UE,Chan,Chan.LSpeed/fc)*...
            ArrayResponse(BS.angle,BS,Chan,Chan.LSpeed/fc)';
        %%% NLoS paths
        for j=2:L
            H_fc(:,:,k)=H_fc(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/fc)*...
            ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/fc)';
        end        
    end
else
    for k=1:N
        f=fc+fs/(N)*(k-1-(N-1)/2);
        %%% NLoS paths
        for j=1:L
            H_fc(:,:,k)=H_fc(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/fc)*...
            ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/fc)';
        end        
    end
end

H_f = H_f.*sqrt(Nt*Nr/L);
H_fc = H_fc.*sqrt(Nt*Nr/L);
% H_fc = H_f;
end