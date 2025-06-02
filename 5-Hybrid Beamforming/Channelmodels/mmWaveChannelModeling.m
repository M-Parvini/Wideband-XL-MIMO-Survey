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
% Number of clusters
Nc = Chan.numClusters;
% Number of rays
Nray = Chan.Rays;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H
tmax = Chan.delay_spread;

%%% Frequency-domain channel
H_f=zeros(Nr,Nt,N);
H_fc=zeros(Nr,Nt,N);

%%% BS and UE AoD and AoA
Alpha_BS = zeros(Nc*Nray, 2);
Alpha_UE = zeros(Nc*Nray, 2);
gains = zeros(Nc*Nray, 1);
delays = rand(Nc*Nray,1)*tmax;

%%% Array responses
Atf  = zeros(Nt, Nc*Nray, N);
Atfc = zeros(Nt, Nc*Nray);
Arf  = zeros(Nr, Nc*Nray, N);
Arfc = zeros(Nr, Nc*Nray);

%%% Channel modeling

for c = 1:Nc
    AoD_m = unifrnd(0,pi,1,2)-pi/2;
    AoA_m = unifrnd(0,pi,1,2)-pi/2;
    % AoD_m = unifrnd(0,2*pi,1,2);
    % AoA_m = unifrnd(0,2*pi,1,2);
    
    AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
    AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
    AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
    AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
    % AoD(1,:) = ones(1,Nray)*AoD_m(1);
    % AoD(2,:) = ones(1,Nray)*AoD_m(2);
    % AoA(1,:) = ones(1,Nray)*AoA_m(1);
    % AoA(2,:) = ones(1,Nray)*AoA_m(2);
    
    %%% Saving the angles
    Alpha_BS((c-1)*Nray+1:c*Nray,:) = AoD.';
    Alpha_UE((c-1)*Nray+1:c*Nray,:) = AoA.';

    for j = 1:Nray
        temp = (c-1)*Nray+j;

        %%% Frequency independent array responses
        Atfc(:,temp) = ArrayResponse([AoD(1,j),AoD(2,j)],BS,Chan,Chan.LSpeed/fc);
        Arfc(:,temp) = ArrayResponse([AoA(1,j),AoA(2,j)],UE,Chan,Chan.LSpeed/fc);
        %%% channel gain
        alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
        % alpha=1;
        gains(temp, 1) = alpha;
        for k=1:N
            fk = fs/(N)*(k-1-(N-1)/2);
            f = fc+fk;

            %%% Frequency dependent array responses
            Atf(:,temp,k) = ArrayResponse([AoD(1,j),AoD(2,j)],BS,Chan,Chan.LSpeed/f);
            Arf(:,temp,k) = ArrayResponse([AoA(1,j),AoA(2,j)],UE,Chan,Chan.LSpeed/f);

            % H_f(:,:,k) = H_f(:,:,k)+alpha*Arf(:,temp)*Atf(:,temp)'*...
            % exp(-1j*2*pi/N*(k-1)*(c-1));
            H_f(:,:,k) = H_f(:,:,k)+alpha*Arf(:,temp,k)*Atf(:,temp,k)'*...
            exp(-1j*2*pi*delays(temp,1)*f);

            H_fc(:,:,k) = H_fc(:,:,k)+alpha*Arfc(:,temp)*Atfc(:,temp)'*...
            exp(-1j*2*pi*delays(temp,1)*f);
            % H_fc(:,:,k) = H_fc(:,:,k)+alpha*Arfc(:,temp)*Atfc(:,temp)'*...
            % exp(-1j*2*pi/N*(k-1)*(c-1));
        end
    end
end

H_f = H_f.*gamma;
H_fc = H_fc.*gamma;

%%% Saving the array response
Chan.Atf  = Atf;
Chan.Atfc = Atfc;
Chan.Arf  = Arf;
Chan.Arfc = Arfc;
Chan.normalizedPathGains = gains;

%%% Saving the BS and UE angles
Chan.BSangle = Alpha_BS(:,1);
Chan.UEangle = Alpha_UE(:,1);

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