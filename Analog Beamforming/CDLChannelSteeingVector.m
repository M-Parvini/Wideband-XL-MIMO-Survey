function [wT, wR] = CDLChannelSteeingVector(Chan, BS, UE)

%%% Steering vector modeling for BS
wT = ArrayResponse(BS.angle,BS,Chan,Chan.lambda);

%%% Steering vector modeling for UE
wR = ArrayResponse(UE.angle,UE,Chan,Chan.lambda);
